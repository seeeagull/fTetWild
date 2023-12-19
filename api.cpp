// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>

#include <Eigen/Dense>
#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>

#ifdef LIBIGL_WITH_TETGEN
#include <igl/copyleft/tetgen/tetrahedralize.h>
#endif

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <bitset>
#include <vector>
#include <cstdlib>

using namespace floatTetWild;
using namespace Eigen;

class GeoLoggerForward : public GEO::LoggerClient
{
    std::shared_ptr<spdlog::logger> logger_;

  public:
    template<typename T>
    GeoLoggerForward(T logger)
        : logger_(logger)
    {}

  private:
    std::string truncate(const std::string& msg)
    {
        static size_t prefix_len = GEO::CmdLine::ui_feature(" ", false).size();
        return msg.substr(prefix_len, msg.size() - 1 - prefix_len);
    }

  protected:
    void div(const std::string& title) override
    {
        logger_->trace(title.substr(0, title.size() - 1));
    }

    void out(const std::string& str) override { logger_->info(truncate(str)); }

    void warn(const std::string& str) override { logger_->warn(truncate(str)); }

    void err(const std::string& str) override { logger_->error(truncate(str)); }

    void status(const std::string& str) override
    {
        // Errors and warnings are also dispatched as status by geogram, but without
        // the "feature" header. We thus forward them as trace, to avoid duplicated
        // logger info...
        logger_->trace(str.substr(0, str.size() - 1));
    }
};

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>
#include <floattetwild/Predicates.hpp>

#include <floattetwild/MshLoader.h>
#include <geogram/mesh/mesh_AABB.h>

inline int str_to_int(const std::string &str) {
    char *end_ptr = nullptr;
    int result = (int) strtol(str.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        throw std::runtime_error("Could not parse signed integer \"" + str + "\"");
    return result;
}

inline float str_to_float(const std::string &str) {
    char *end_ptr = nullptr;
    float result = (float) strtod(str.c_str(), &end_ptr);
    if (*end_ptr != '\0')
        throw std::runtime_error("Could not parse floating point value \"" + str + "\"");
    return result;
}

int runFTetWild(std::vector<std::vector<int>> &faces,
                std::vector<std::vector<float>> &verts,
                int argc, char **argv) {
#ifdef STORE_SAMPLE_POINTS
    cout << "STORE_SAMPLE_POINTS defined" << endl;
#endif

#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    GEO::initialize();
    //    exactinit();

    // Import standard command line arguments, and custom ones
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    bool run_tet_gen   = false;
    bool skip_simplify = false;
    bool nobinary      = false;
    bool nocolor       = false;
    bool export_raw    = false;
    int boolean_op = -1;
    std::string csg_file   = "";
    std::string background_mesh = "";
    std::string epsr_tags = "";

    Mesh        mesh;
    Parameters& params = mesh.params;

    try {
        for (int i=1; i<argc; ++i) {
            if (strcmp("-i", argv[i]) == 0) {
                params.input_path = argv[++i];
            } else if (strcmp("-o", argv[i]) == 0) {
                params.output_path = argv[++i];
            } else if (strcmp("--tag", argv[i]) == 0) {
                params.tag_path = argv[++i];
            } else if (strcmp("--op", argv[i]) == 0) {
                boolean_op = str_to_int(argv[++i]);
            } else if (strcmp("-e", argv[i]) == 0) {
                params.eps_rel = str_to_float(argv[++i]);
            } else if (strcmp("-l", argv[i]) == 0) {
                params.ideal_edge_length_rel = str_to_float(argv[++i]);
            } else if (strcmp("--stop-energy", argv[i]) == 0) {
                params.stop_energy = str_to_float(argv[++i]);
            } else if (strcmp("--skip-simplify", argv[i]) == 0) {
                skip_simplify = true;
            } else if (strcmp("--no-binary", argv[i]) == 0) {
                nobinary = false;
            } else if (strcmp("--no-color", argv[i]) == 0) {
                nocolor = true;
            } else if (strcmp("--export-raw", argv[i]) == 0) {
                export_raw = true;
            } else if (strcmp("--smooth-open-boundary", argv[i]) == 0) {
                params.smooth_open_boundary = true;
            } else if (strcmp("--manifold-surface", argv[i]) == 0) {
                params.manifold_surface = true;
            } else if (strcmp("--coarsen", argv[i]) == 0) {
                params.coarsen = true;
            } else if (strcmp("--csg", argv[i]) == 0) {
                csg_file = argv[++i];
            } else if (strcmp("--disable-filtering", argv[i]) == 0) {
                params.disable_filtering = true;
            } else if (strcmp("--use-floodfill", argv[i]) == 0) {
                params.use_floodfill = true;
            } else if (strcmp("--use-general-wn", argv[i]) == 0) {
                params.use_general_wn = true;
            } else if (strcmp("--use-input-for-wn", argv[i]) == 0) {
                params.use_input_for_wn = true;
            } else if (strcmp("--bg-mesh", argv[i]) == 0) {
                background_mesh = argv[++i];
            } else if (strcmp("--level", argv[i]) == 0) {
                params.log_level = str_to_int(argv[++i]);
            }
        }
    } catch (const std::exception &e) {
        cout << "Error: " << e.what() << endl;
    }

#ifdef LIBIGL_WITH_TETGEN
    command_line.add_flag("--tetgen", run_tet_gen, "run tetgen too. (optional)");
#endif
    unsigned int max_threads = std::numeric_limits<unsigned int>::max();

#ifdef FLOAT_TETWILD_USE_TBB
    const size_t MB          = 1024 * 1024;
    const size_t stack_size  = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(max_threads, num_threads);
    params.num_threads       = num_threads;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

    Logger::init(!params.is_quiet, params.log_path);
    params.log_level = std::max(0, std::min(6, params.log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    GEO::Logger* geo_logger = GEO::Logger::instance();
    geo_logger->unregister_all_clients();
    geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
    geo_logger->set_pretty(false);

    if (!verts.empty())
        params.input_path.clear();
    if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    std::string output_mesh_name = params.output_path;
    if (params.output_path.size() > 3 &&
        params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) ==
          "msh")
        output_mesh_name = params.output_path;
    else if (params.output_path.size() > 4 &&
             params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) ==
               "mesh")
        output_mesh_name = params.output_path;
    else if (params.output_path.size() > 0)
        output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// set sizing field
    Eigen::VectorXd V_in;
    Eigen::VectorXi T_in;
    Eigen::VectorXd values;
    if (!background_mesh.empty()) {
        PyMesh::MshLoader mshLoader(background_mesh);
        V_in   = mshLoader.get_nodes();
        T_in   = mshLoader.get_elements();
        values = mshLoader.get_node_field("values");
    }
    if (V_in.rows() != 0 && T_in.rows() != 0 && values.rows() != 0) {
        params.apply_sizing_field = true;

        params.V_sizing_field = V_in;
        params.T_sizing_field = T_in;
        params.values_sizing_field = values;
    }

    /// set input tage
    std::vector<Eigen::Vector3d>  input_vertices{};
    std::vector<Eigen::Vector3i> input_faces{};
    std::vector<int>      input_tags{};
    //    std::vector<double> input_epsr_tags{};

    if (!params.tag_path.empty()) {
        input_tags.reserve(input_faces.size());
        std::string   line;
        std::ifstream fin(params.tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                input_tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty()) {
        std::ifstream fin(epsr_tags);
        std::string   line;
        while (std::getline(fin, line)) {
            params.input_epsr_tags.push_back(std::stod(line));
        }
        fin.close();
    }
#endif

    /// set envelope
    igl::Timer               timer;
    GEO::Mesh                sf_mesh;
    json                     tree_with_ids;
    std::vector<std::string> meshes;
    if (!csg_file.empty()) {
        json          csg_tree = json({});
        std::ifstream file(csg_file);

        if (file.is_open())
            file >> csg_tree;
        else {
            logger().error("unable to open {} file", csg_file);
            return EXIT_FAILURE;
        }
        file.close();

        CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

        if (!CSGTreeParser::load_and_merge(
              meshes, input_vertices, input_faces, sf_mesh, input_tags))
            return EXIT_FAILURE;

        // To disable the recent modification of using input for wn, use meshes.clear();
    }
    else {
        if (!verts.empty()) {
            if (!MeshIO::construct_mesh(verts, faces, input_vertices, input_faces, sf_mesh, input_tags)) {
                return EXIT_FAILURE;
            }
        } else {
#ifdef NEW_ENVELOPE
            if (!MeshIO::load_mesh(params.input_path,
                                input_vertices,
                                input_faces,
                                sf_mesh,
                                input_tags,
                                params.input_epsr_tags)) {
#else
            if (!MeshIO::load_mesh(
                params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
#endif
                logger().error("Unable to load mesh at {}", params.input_path);
                MeshIO::write_mesh(output_mesh_name, mesh, false);
                return EXIT_FAILURE;
            }
            else if (input_vertices.empty() || input_faces.empty()) {
                MeshIO::write_mesh(output_mesh_name, mesh, false);
                return EXIT_FAILURE;
            }
        }

        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
    }
    AABBWrapper tree(sf_mesh);
    if (!params.init(tree.get_sf_diag())) {
        return EXIT_FAILURE;
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty())
        tree.init_sf_tree(
          input_vertices, input_faces, params.input_epsr_tags, params.bbox_diag_length);
    else
        tree.init_sf_tree(input_vertices, input_faces, params.eps);
#endif

#ifdef LIBIGL_WITH_TETGEN
    if (run_tet_gen) {
        Eigen::MatrixXd tetgen_pts(input_vertices.size(), 3);
        Eigen::MatrixXi tetgen_faces(input_faces.size(), 3);

        for (size_t i = 0; i < input_vertices.size(); ++i) {
            tetgen_pts.row(i) = input_vertices[i].cast<double>();
        }

        for (size_t i = 0; i < input_faces.size(); ++i) {
            tetgen_faces.row(i) = input_faces[i];
        }

        std::stringstream buf;
        buf.precision(100);
        buf.setf(std::ios::fixed, std::ios::floatfield);
        buf << "Qpq2.0a"
            << params.ideal_edge_length * params.ideal_edge_length * params.ideal_edge_length *
                 sqrt(2.) / 12.;

        Eigen::MatrixXi tetgen_generated_tets;
        Eigen::MatrixXd tetgen_generated_points;
        Eigen::MatrixXi tetgen_generated_faces;

        timer.start();
        igl::copyleft::tetgen::tetrahedralize(tetgen_pts,
                                              tetgen_faces,
                                              buf.str(),
                                              tetgen_generated_points,
                                              tetgen_generated_tets,
                                              tetgen_generated_faces);
        timer.stop();
        logger().info("Tetgen time {}s", timer.getElapsedTimeInSec());
        stats().record(StateInfo::tetgen_id,
                       timer.getElapsedTimeInSec(),
                       tetgen_generated_points.rows(),
                       tetgen_generated_tets.rows(),
                       0,
                       0);
    }
#endif

    stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

    timer.start();
    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
    logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::preprocessing_id,
                   timer.getElapsedTimeInSec(),
                   input_vertices.size(),
                   input_faces.size(),
                   -1,
                   -1);
    if (params.log_level <= 1)
        output_component(input_vertices, input_faces, input_tags);

    timer.start();
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::tetrahedralization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   -1,
                   -1);

    timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    logger().info("cutting {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::cutting_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy(),
                   std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    timer.start();
    optimization(
      input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
    logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::optimization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());

    timer.start();
    correct_tracked_surface_orientation(mesh, tree);
    logger().info("correct_tracked_surface_orientation done");

    if (export_raw) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 3> Vt;
        Eigen::Matrix<int, Eigen::Dynamic, 3>    Ft;

        if (!csg_file.empty()) {
            int max_id = CSGTreeParser::get_max_id(tree_with_ids);

            for (int i = 0; i <= max_id; ++i) {
                get_tracked_surface(mesh, Vt, Ft, i);
                igl::write_triangle_mesh(
                  params.output_path + "_" + params.postfix + "_" + std::to_string(i) + "_all.obj",
                  Vt,
                  Ft);
            }
        }
        else {
            get_tracked_surface(mesh, Vt, Ft);
            igl::write_triangle_mesh(
              params.output_path + "_" + params.postfix + "_all.obj", Vt, Ft);
        }
        MeshIO::write_mesh(params.output_path + "_" + params.postfix + "_all.msh", mesh, false);
    }

    if (!csg_file.empty())
        boolean_operation(mesh, tree_with_ids, meshes);
    else if (boolean_op >= 0)
        boolean_operation(mesh, boolean_op);
    else {
        if (params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto& t : mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        }
        else {
            if (!params.disable_filtering) {
                if (params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                }
                else if (params.use_input_for_wn) {
                    filter_outside(mesh, input_vertices, input_faces);
                }
                else
                    filter_outside(mesh);
            }
        }
    }
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    if (params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    stats().record(StateInfo::wn_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());
    logger().info("after winding number");
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("winding number {}s", timer.getElapsedTimeInSec());
    logger().info("");

    // fortest
    std::vector<Scalar> colors;
    if (!nocolor) {
        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
    }
    // fortest
    if (params.output_path.size() > 0) {
        MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !nobinary, !csg_file.empty());
        igl::write_triangle_mesh(params.output_path + "_" + params.postfix + "_sf.obj", V_sf, F_sf);
    }
    MeshIO::extract_mesh_element(mesh, faces, verts);

    if (!params.log_path.empty()) {
        std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
        if (fout.good())
            fout << stats();
        fout.close();
        if (!params.envelope_log.empty()) {
            std::ofstream fout(params.envelope_log);
            fout << envelope_log_csv;
            fout.close();
        }
    }

    return EXIT_SUCCESS;
}