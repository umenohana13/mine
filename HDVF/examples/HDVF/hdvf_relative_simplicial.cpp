#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_2.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Hdvf_relative.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_2<Kernel> Traits;

typedef Kernel::Point_2 Point_2;

//typedef CGAL::Zp<5,int,true> Coefficient_ring;
typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sub_sparse_matrix<CGAL::OSM::Sparse_chain> Sparse_matrix_struct;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits, Sparse_matrix_struct> Complex;
typedef HDVF::Hdvf_relative<Complex> HDVF_type;
typedef HDVF::Sub_chain_complex_mask<Complex> Sub_chain_complex;

int main(int argc, char **argv)
{
    std::string filename, subset_root_name, nodes_file ;
    size_t n_subs;
    if (argc != 5) std::cout << "usage: hdvf_relative_simplicial simp_file subset_root_name number_of_subsets nodes_file" << std::endl;
    else {
        filename = argv[1];
        subset_root_name = argv[2];
        n_subs = std::stoi(argv[3]);
        nodes_file = argv[4];

        // Load mesh object
        HDVF::Mesh_object_io<Traits> mesh ;
        mesh.read_simp(filename);
        mesh.read_nodes_file(nodes_file);
        mesh.print_infos();

        // Build the complex
        Complex complex(mesh);

        std::cout << "--------- Flow cell 1 ---------" << std::endl;

        // Load first Sub_chain_complex_mask;
        std::string sub_file_name(subset_root_name+"1.sub");
        HDVF::Sub_chain_complex_mask<Complex> K_init(CGAL::IO::read_SUB(complex, sub_file_name, false));

        std::cout << "------> subcomplex K" << std::endl ;
        std::cout << K_init << std::endl ;

        // Create and compute a perfect HDVF over the complex with mask K_init
        HDVF_type hdvf(complex, K_init, HDVF::OPT_BND, false);
        hdvf.compute_perfect_hdvf();

        CGAL::IO::write_VTK(hdvf, complex, "tmp/subset1", false);

        std::vector<std::vector<size_t> > criticals(hdvf.psc_flags(HDVF::CRITICAL));
        for (int q=0; q<= complex.dimension(); ++q) {
            std::cout << "criticals dim " << q << ": ";
            for (size_t i : criticals.at(q))
                std::cout << i << " ";
            std::cout << std::endl;
        }

        // Compute perfect HDVFs over next subsets
        for (int i = 2; i < n_subs; ++i) {
            std::cout << "--------- Flow cell " << std::to_string(i) << " ---------" << std::endl;
            std::string sub_file_name(subset_root_name+std::to_string(i)+".sub");
            HDVF::Sub_chain_complex_mask<Complex> K(CGAL::IO::read_SUB(complex, sub_file_name, false));
            hdvf.change_subset_complex(K, false);
            hdvf.compute_perfect_hdvf();

            std::vector<std::vector<size_t> > criticals(hdvf.psc_flags(HDVF::CRITICAL));
            for (int q=0; q<= complex.dimension(); ++q) {
                std::cout << "criticals dim " << q << ": ";
                for (size_t i : criticals.at(q))
                    std::cout << i << " ";
                std::cout << std::endl;
            }

            CGAL::IO::write_VTK(hdvf, complex, "tmp/subset"+std::to_string(i), false);
        }

        // Full (non perfect) HDVF
        HDVF::Sub_chain_complex_mask<Complex> K_full(complex);
        hdvf.change_subset_complex(K_full,true, false);
        CGAL::IO::write_VTK(hdvf, complex, "tmp/full", false);

        hdvf.print_reduced_boundary();
    }

        return 0;
}
