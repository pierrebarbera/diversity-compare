/*
    Copyright (C) 2019 Pierre Barbera

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Pierre Barbera <pierre.barbera@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <utility>
#include <tuple>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <vector>
#include <cmath>
#include <math.h>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

bool equals_approx( double const a, double const b , double const epsilon=1e-10)
{
    return std::abs( a - b ) < epsilon;
}

double phylo_entropy_g( double const x )
{
    if( equals_approx( x, 0.0 ) ) {
        return 0.0;
    }
    return x * std::log(x);
}

double phylo_quad_entropy_g( double const x )
{
    return x * (1.0 - x);
}

double step_function_g( double const x, double const theta )
{
    assert( equals_approx(theta, 1.0) or equals_approx(theta, 0.0) or (theta > 0.0 and theta < 1.0) );
    assert( equals_approx(x, 1.0) or equals_approx(x, 0.0) or (x > 0.0 and x < 1.0) );

    // special case: only consider edges in the spanning tree of the sample
    // (0 or 1 fraction means non-shared edge, which means the edge is not in the spanning)
    if( equals_approx( x, 0.0 ) or equals_approx( x, 1.0 ) ) {
        return 0.0;
    }

    return std::pow( 2 * std::min( x, 1.0 - x ), theta );
}

size_t num_queries( std::vector< Pquery const* > const& pqs )
{
    size_t sum = 0;
    for( auto const& pq_ptr : pqs ) {
        sum += pq_ptr->name_size();
    }
    return sum;
}

size_t num_queries( Sample const& sample )
{
    return total_name_count( sample );
}

double mass_prox_dist( MassTreeEdge const& edge )
{
    auto const& masses = edge.data< MassTreeEdgeData >().masses;
    double most_distant = 0.0;
    for (auto it = masses.begin(); it != masses.end(); ++it) {
        most_distant = ( most_distant < it->first ) ? it->first : most_distant;
    }
    return most_distant;
}

/**
 * BWPD-style branch length sum of mass-edges (made by placements for example)
 * on the given edge of the reference tree.
 */
template< class g_func_t, class ...Args >
double distal_edge_sum( MassTreeEdge const& edge, g_func_t g_func, Args... args )
{
    double sum = 0.0;
    double dragged_mass = 0.0;
    auto const& masses = edge.data< MassTreeEdgeData >().masses;

    // iterate in reverse order through the map of masses, starting with the most distal mass
    for( auto it = masses.rbegin(); it != masses.rend(); ++it ) {
        auto prox_length    = it->first;

        // the "distal" mass as it gets dragged along the path of masses (a monofurcating subtree of sorts)
        dragged_mass += it->second;

        assert( dragged_mass > 0.0 or equals_approx( dragged_mass, 0.0 ) );
        assert( dragged_mass < 1.0 or equals_approx( dragged_mass, 1.0 ) );

        double branch_length = prox_length;

        // look ahead (or rather, behind from the perspective of the container) to get the next
        // proximal length, as <current prox length> - <next prox length> = "branch length"
        // of the mass
        auto next_it = std::next( it );
        if( next_it != masses.rend() ) {
            auto next_prox_length = next_it->first;
            branch_length = prox_length - next_prox_length;
        }

        assert( branch_length > 0.0 );

        // add the sum according to the function
        sum += branch_length * g_func( dragged_mass, args... );

    }
    return sum;
}

template< class g_func_t, class ...Args >
double MassTreeBWPD( MassTree const& mass_tree, g_func_t g_func, Args... args )
{
    /*  The goal is to iterate over all edges and calculate D_s(i), which is
        the fraction of reads in sample s that are in leaves (or edges / placements
        in our case) on the distal side of edge i. Then, the BWPD is the sum of edge
        length l(i) multiplied by [2min(D(i),1−D(i))]^θ. For θ=0, and evaluating only
        those edges that are in a samples spanning tree, this results in the somewhat
        classical FaithPD, except adapted to phylogenetic placement.
    */

    // we start by first calculating D(i) for every edge in the tree
    // which we do bottom-up, post-order, dragging with us previous distal-side
    // counts, dividing them by the total count

    auto mass_per_edge = mass_tree_mass_per_edge( mass_tree );

    std::vector< double > per_edge_D( mass_tree.edge_count(), -1.0 );

    double result = 0.0;

    size_t count = 0;

    // traverse the tree
    for( auto const& it : postorder( mass_tree ) ) {
        // ensure last edge isn't visited twice
        if( it.is_last_iteration() ) { continue; }
        ++count;

        auto const& edge = it.edge();
        auto const& edge_index = edge.index();

        // if this is a leaf, there cannot be any mass on the distal side, so we skip
        if( is_leaf( edge ) ) {
            per_edge_D[ edge_index ] = 0.0;
            continue;
        }
        // interior edges:
        // calculate the new distal count as the sum of distal counts and query counts
        // of the child edges
        auto const& node = it.node();
        auto const& lhs_edge    = node.link().next().edge();
        auto const& rhs_edge    = node.link().next().next().edge();
        auto const lhs_edge_index = lhs_edge.index();
        auto const rhs_edge_index = rhs_edge.index();

        per_edge_D[ edge_index ] =
            per_edge_D[ lhs_edge_index ] + mass_per_edge[ lhs_edge_index ] +
            per_edge_D[ rhs_edge_index ] + mass_per_edge[ rhs_edge_index ];

        auto branch_length = edge.data< MassTreeEdgeData >().branch_length;

        // "recurse" to get the part of the sum coming from the basses located on the immediate
        // child edges, which are treated differently since they are not directly in the tree
        // structure
        if( equals_approx( per_edge_D[ lhs_edge_index ], 0.0 ) ){
            result += distal_edge_sum( lhs_edge, g_func, std::forward<Args>(args)... );
        }
        if( equals_approx( per_edge_D[ rhs_edge_index ], 0.0 ) ){
            result += distal_edge_sum( rhs_edge, g_func, std::forward<Args>(args)... );
        }

        // update the BWPD sum with the result of this edge
        result += branch_length * g_func( per_edge_D[ edge_index ], args... );
    }

    if( count != mass_tree.edge_count() ) {
        throw std::runtime_error("count != mass_tree.edge_count()");
    }

    return result;
}

MassTree convert_key_attribute_tree_to_scrapp_mass_tree( AttributeTree const& source )
{
    return convert(
        source,
        [] ( tree::BaseNodeData const& node_data ) {
            auto new_node_data = tree::MassTreeNodeData::create();
            auto const& orig_node = dynamic_cast< AttributeTreeNodeData const& >( node_data );
            new_node_data->name = orig_node.name;

            return new_node_data;
        },
        [] ( tree::BaseEdgeData const& edge_data ) {
            auto new_edge_data = tree::MassTreeEdgeData::create();
            auto const& orig_edge = dynamic_cast< AttributeTreeEdgeData const& >( edge_data );
            new_edge_data->branch_length = orig_edge.branch_length;
            auto bl_half = orig_edge.branch_length / 2.0;
            auto weight = stod( orig_edge.attributes.at( "species_count" ) );
            assert( orig_edge.attributes.count( "species_count" ) > 0 );
            new_edge_data->masses[ bl_half ] = weight;

            return new_edge_data;
        }
    );
}

/**
 *  Calculate different diversity metrics based on the given jplace files, output as CSV
 */
int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if ( argc < 2 ) {
        throw std::runtime_error(
            std::string("Usage: ") + argv[0] + " <scrapp-files...>\n"
        );
    }

    std::vector<std::string> scrapp_files;
    for ( int i = 1; i < argc; ++i ) {
        scrapp_files.emplace_back(argv[i]);
    }

    std::vector< double > theta_set = { 0.0, 0.25, 0.5, 0.75, 1.0 };

    // write the header
    std::cout << "sample,phylo_entropy,quadratic,bwpd";
    for( auto const theta : theta_set ) {
        std::cout << ",bwpd_" << theta;
    }
    std::cout << "\n";

    // set up the reader to parse NHX attributes, specifically the ones we want
    auto reader = KeyedAttributeTreeNewickReader();
    reader.set_nhx_delimiters();
    reader.add_attribute( "species_count", KeyedAttributeTreeNewickReaderPlugin::Target::kEdge, "species_count", "0.0" );

    for( size_t i = 0; i < scrapp_files.size(); ++i ) {
        auto attr_tree = reader.read( from_file( scrapp_files[ i ] ) );
        auto mass_tree = convert_key_attribute_tree_to_scrapp_mass_tree( attr_tree );
        mass_tree_normalize_masses( mass_tree );

        std::cout << file_basename( file_path( scrapp_files[ i ] ) );

        if( not is_bifurcating( mass_tree ) ) {
            throw std::runtime_error("non bifurcating input tree!");
        }

        // Phylogenetic Entropy
        std::cout  << "," << -MassTreeBWPD( mass_tree, phylo_entropy_g );

        // Phylogenetic Quadratic Entropy
        std::cout  << "," << MassTreeBWPD( mass_tree, phylo_quad_entropy_g );

        // BWPD
        std::cout  << "," << MassTreeBWPD( mass_tree, step_function_g, 1.0 );
        for( auto const theta : theta_set ) {
            std::cout  << "," << MassTreeBWPD( mass_tree, step_function_g, theta );
        }
        std::cout << "\n";
    }

    return 0;
}
