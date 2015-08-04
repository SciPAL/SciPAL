#ifndef DRS_CELL_HH
#define DRS_CELL_HH

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

// Given the complexity of the dielectric relaxation problem
// we put everything into a dedicated namespace.
namespace DRS {

// @sect3{struct BCids}
//
// The interesting things of the DRS problem happen on the boundaries.
// In order to address them by name we define an auxiliary structure
// which provides the mapping from names to numbers such that
// we can use deal.II's subboundary facilites.
struct BCids {

    static unsigned char anode()            { return 5; }

    static unsigned char cylinder_surface() { return 0; }

    static unsigned char cavity()           { return 3; }

    static unsigned char cathode()          { return 1; }
};

// @sect3{class MeshGenerator}
//
// The coarse mesh needs quite a few parameters. We pack the parameters into
// a structure @p DRSCellParams which manages
// setting the parameters at run-time using deal.II's ParameterHandler class.
// The struct for the parameters and the function @p drs_cell which actually
// generates the mesh are bundled in a class @p MeshGenerator.
class MeshGenerator {

public:
    struct DRSCellParams {

        static void declare(dealii::ParameterHandler & prm)
        {
            prm.enter_subsection("DRS geometry");

            prm.declare_entry("Cavity diameter", "1.",
                              dealii::Patterns::Double(),
                              "Radius of dielectric sample in units of Debye length.");

            prm.declare_entry("Cylinder diameter", "10.",
                              dealii::Patterns::Double(),
                              "In units of Debye length.");

            prm.declare_entry("Cylinder height", "10.",
                              dealii::Patterns::Double(),
                              "In units of Debye length.");

            prm.declare_entry("Print VX", "false",
                              dealii::Patterns::Bool() ,
                              "For debugging the geometry it is helpful to be able to print "
                              "all the vertices of the coarse mesh.");

            prm.leave_subsection();
        }

        void get(dealii::ParameterHandler & prm)
        {
            prm.enter_subsection("DRS geometry");

            cavity_radius = prm.get_double("Cavity diameter")/2.;

            cyl_radius = prm.get_double("Cylinder diameter")/2.;

            cyl_half_height = prm.get_double("Cylinder height")/2.;

            print_vx = prm.get_bool("Print VX");

            prm.leave_subsection();
        }

        double cavity_radius, cyl_radius, cyl_half_height;

        bool print_vx;
    };

    template<int dim>
    static void drs_cell (dealii::Triangulation<dim>& tria, const DRSCellParams &params, //std::string prm_drs_cell_filepath,
                          int n_init_refinements)
    {
        typedef  dealii::Point<dim> Point;



        // Once the parameters are read we log what we actually have read.

        //            MakeParams<DRSCellParams>::initialize_prm_handler(prm,
        //            (PathInfo::prmfile_homedir / "drs_cell_geo").string(),
        //            (PathInfo::rundir /
        //            "log").string() );



        {
            // Implementation for 3D only
            Assert (dim == 3, dealii::ExcNotImplemented());

            /*
                ----------------------------------------
               |            |                |          |
               |----------------------------------------|
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |----------------------------------------|
               |            \               /           |
               |             \             /            |
               |              \           /             |
               |-----          \         /         -----|
               |      -----     \       /     -----     |
               |           -----  -----  -----          |
               |                 |     |                |
               |                 |     |                |
               |           -----  -----  -----          |
               |     -----      /       \     -----     |
               |----           /         \         -----|
               |              /           \             |
               |             /             \            |
               |            /               \           |
               |----------------------------------------|
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |            |                |          |
               |----------------------------------------|
               |            |                |          |
                ----------------------------------------

              -d           -a               +a         +d
    */

            // Copy the base from hyper_ball<dim>
            // and transform it to yz
            const double d = params.cyl_radius/std::sqrt(2.0);
            const double a = d/(1+std::sqrt(2.0));

            // To create the coarse mesh we need a few magic numbers.
#define n_vertices_per_cell  8
#define n_vertices_per_layer 8
#define n_cells_per_layer    5
#define n_cell_layers        7

            static const unsigned int n_vertices = n_vertices_per_layer * (n_cell_layers + 1);

            const unsigned int cavity_layer_id = n_cell_layers/2;

            typedef dealii::Point<dim> Point;

            Point layer_vertices[n_vertices_per_layer] = {
                Point(-d, 0,-d),
                Point( d, 0,-d),
                Point(-a, 0,-a),
                Point( a, 0,-a),
                Point(-a, 0, a),
                Point( a, 0, a),
                Point(-d, 0, d),
                Point( d, 0, d)
            };

            double height[n_cell_layers+1] = {
                0.,
                1.,
                5.,
                8.,
                12.,
                15.,
                19.,
                20.
            };

            // Rescale to desired height
            for (unsigned int i = 0; i < n_cell_layers+1; i++)
                height[i] *= 2*params.cyl_half_height/20.;

            const double a_cav = params.cavity_radius/std::sqrt(3.);
            const double cav_center = 0.;
            Point cavity_layer_shifts[4] = {
                Point(-a_cav, cav_center, -a_cav),
                Point( a_cav, cav_center, -a_cav),
                Point(-a_cav, cav_center,  a_cav),
                Point( a_cav, cav_center,  a_cav)
            };

            const double max_height = height[n_cell_layers];

            std::vector<Point> vertices(n_vertices);

            typename std::vector<Point>::iterator v = vertices.begin(); // global vertex counter;

            for (unsigned int l = 0; l < n_cell_layers+1; l++)
            {
                for (unsigned int vl = 0; vl < n_vertices_per_layer; vl++)
                {
                    Point & vertex = *v;
                    vertex    = layer_vertices[vl];
                    vertex(1) = height[l]-max_height/2;

                    // shift cavity corners
                    if ((l == cavity_layer_id) && ( vl >=2) && (vl < 6) )
                    {
                        vertex = cavity_layer_shifts[vl-2];
                        vertex(1) -= a_cav;
                    }

                    if ((l == cavity_layer_id+1) && ( vl >=2) && (vl < 6) )
                    {
                        vertex = cavity_layer_shifts[vl-2];
                        vertex(1) += a_cav;
                    }

                    // Turn cylinder such that y becomes x
                    {
                        const double h = vertex(1);
                        vertex(1) = -vertex(0);
                        vertex(0) = h;
                    }

                    if (params.print_vx )
                        std::cout << vertex << std::endl;
                    ++v;
                }
            }

            int layer_cell_vertices[n_cells_per_layer][n_vertices_per_cell] = {
                {0, 1, 8, 9, 2, 3, 10, 11},
                {0, 2, 8, 10, 6, 4, 14, 12},
                {2, 3, 10, 11, 4, 5, 12, 13},
                {1, 7, 9, 15, 3, 5, 11, 13},
                {6, 4, 14, 12, 7, 5, 15, 13}
            };

            // The cavity is a non-existing cell, thus we have to subtract one.
            const unsigned int n_cells = n_cells_per_layer * n_cell_layers-1;

            std::vector<dealii::CellData<dim> > cells (n_cells, dealii::CellData<dim>());

            // We need a global cell counter.
            unsigned int c = 0;
            // Then, we loop over the layers from bottom to top.
            for (unsigned int l = 0; l < n_cell_layers; l++)
                for (unsigned int cl=0; cl<n_cells_per_layer; ++cl)
                {
                    if ((l == cavity_layer_id) && (cl == 2))
                        continue;
                    else {
                        for (unsigned int j=0; j<n_vertices_per_cell; ++j)
                            cells[c].vertices[j] = layer_cell_vertices[cl][j]+l*n_vertices_per_layer;
                        cells[c].material_id = 0;
                        c++;
                    }
                };

            // At the beginning, we create a triangulation without any boundary information.
            tria.create_triangulation (vertices, cells,
                                       dealii::SubCellData());

            // The boundary indicators for the
            // faces at the ends are set to 1 and 2,
            // respectively. Note that we also
            // have to deal with those lines
            // that are purely in the interior
            // of the ends. We determine wether
            // an edge is purely in the
            // interior if one of its vertices
            // is at coordinates '+-a' as set
            // above.
            typename dealii::Triangulation<dim>::cell_iterator
                    cell =  tria.begin(),
                    end =  tria.end();

            const double half_length = max_height/2;
            for (; cell != end; ++cell)
                for (unsigned int i=0; i<dealii::GeometryInfo<dim>::faces_per_cell; ++i)
                    if (cell->at_boundary(i))
                    {
                        if (cell->face(i)->center()(0) > half_length-1.e-5)
                        {
                            // top
                            cell->face(i)->set_boundary_indicator(BCids::anode() );

                            for (unsigned int e=0; e<dealii::GeometryInfo<dim>::lines_per_face; ++e)
                                if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
                                    cell->face(i)->line(e)->set_boundary_indicator(BCids::anode() );
                        }
                        else if (cell->face(i)->center()(0) < -half_length+1.e-5)
                        {
                            // bottom
                            cell->face(i)->set_boundary_indicator(BCids::cathode() );

                            for (unsigned int e=0; e<dealii::GeometryInfo<dim>::lines_per_face; ++e)
                                if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
                                    cell->face(i)->line(e)->set_boundary_indicator(BCids::cathode() );
                        }
                        else if ( cell->face(i)->center().norm() < 1.001*params.cavity_radius )
                        {
                            // cavity
                            cell->face(i)->set_boundary_indicator(BCids::cavity() );
                            for (unsigned int e=0; e<dealii::GeometryInfo<dim>::lines_per_face; ++e)
                                if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a_cav) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a_cav) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a_cav) ||
                                        (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a_cav))
                                    cell->face(i)->line(e)->set_boundary_indicator(BCids::cavity() );
                        }
                    }

            // The final step is to assign the boundary objects.
            static dealii::CylinderBoundary<dim> cyl(params.cyl_radius);

            static dealii::HyperBallBoundary<dim> cav(dealii::Point<dim>(), params.cavity_radius);

            // cylinder hull
            tria.set_boundary(BCids::cylinder_surface(), cyl);

            tria.set_boundary(BCids::cavity(), cav);

            // Once the coarse mesh is set up we refine it a bit
            // such that the initial computataions is not completely
            // way of from the needed spatial resolution.
            tria.refine_global(n_init_refinements);

            std::ofstream output_file((// PathInfo::rundir /
                                       "drs_cell_coarse_grid.msh") );

            dealii::GridOut().write_msh ( tria, output_file);

            std::cout << "DRS cell is set up." << std::endl;
        }
    }
};

} // END namespace DRS

#endif // DRS_CELL_HH
