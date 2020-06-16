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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

bool equals_approx( double const a, double const b, double const epsilon = 1e-10 )
{
  return std::abs( a - b ) < epsilon;
}

double phylo_entropy_g( double const x )
{
  if( equals_approx( x, 0.0 ) ) {
    return 0.0;
  }
  return x * std::log( x );
}

double phylo_quad_entropy_g( double const x )
{
  return x * ( 1.0 - x );
}

double step_function_g( double const x, double const theta )
{
  assert( theta >= 0.0 and theta <= 1.0 );

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

/**
 * https://dx.doi.org/10.7717%2Fpeerj.157
 *
 * @param  sample [description]
 * @param  theta  [description]
 * @return        [description]
 */
template< class g_func_t, class... Args >
double BWPD( Sample const& sample, g_func_t g_func, Args... args )
{
  /*  The goal is to iterate over all edges and calculate D_s(i), which is
        the fraction of reads in sample s that are in leaves (or edges / placements
        in our case) on the distal side of edge i. Then, the BWPD is the sum of edge
        length l(i) multiplied by [2min(D(i),1−D(i))]^θ. For θ=0, and evaluating only
        those edges that are in a samples spanning tree, this results in the somewhat
        classical FaithPD, except adapted to phylogenetic placement.
    */

  // assert( theta >= 0.0 and theta <= 1.0 );

  // we start by first calculating D(i) for every edge in the tree
  // which we do bottom-up, post-order, dragging with us previous distal-side
  // counts, dividing them by the total count

  auto const& place_tree = sample.tree();
  auto pqrs_per_edge     = pqueries_per_edge( sample, true );

  // if we only take the best hits, num_placements = num_pqueries, making the
  // total size, including multiplicities, simply the name count:
  double const total_queries = num_queries( sample );

  // vector to track the number of queries on the distal side of the edge
  std::vector< size_t > per_edge_distal_count( place_tree.edge_count() );
  // same idea, but for D(i)
  std::vector< double > per_edge_D( place_tree.edge_count() );

  double result = 0.0;

  // traverse the tree
  for( auto const& it : postorder( place_tree ) ) {
    // ensure last edge isn't visited twice
    if( it.is_last_iteration() ) {
      continue;
    }

    auto const& edge       = it.edge();
    auto const& edge_index = edge.index();

    // if this is a leaf, there cannot be any mass on the distal side, so we skip
    if( is_leaf( edge ) ) {
      per_edge_distal_count[ edge_index ] = 0;
      per_edge_D[ edge_index ]            = 0.0;
      // per_edge_distal_count[ edge_index ]   = num_queries( pqrs_per_edge[ edge_index ] );
      // per_edge_D[ edge_index ]              = per_edge_distal_count[ edge_index ] / total_queries;
      continue;
    }

    // interior edges:
    // calculate the new distal count as the sum of distal counts and query counts
    // of the child edges
    auto const& node          = it.node();
    auto const& lhs_edge      = node.link().next().edge();
    auto const& rhs_edge      = node.link().next().next().edge();
    auto const lhs_edge_index = lhs_edge.index();
    auto const rhs_edge_index = rhs_edge.index();

    per_edge_distal_count[ edge_index ] = per_edge_distal_count[ lhs_edge_index ] + num_queries( pqrs_per_edge[ lhs_edge_index ] ) + per_edge_distal_count[ rhs_edge_index ] + num_queries( pqrs_per_edge[ rhs_edge_index ] );

    // per_edge_distal_count[ edge_index ] =
    //     per_edge_distal_count[ lhs_edge_index ] +
    //     per_edge_distal_count[ rhs_edge_index ] + num_queries( pqrs_per_edge[ edge_index ] );

    // calculate D(i) for this edge
    per_edge_D[ edge_index ] = per_edge_distal_count[ edge_index ] / total_queries;

    auto const branch_length = edge.data< PlacementEdgeData >().branch_length;

    // update the BWPD sum with the result of this edge
    result += branch_length * g_func( per_edge_D[ edge_index ], args... );
  }

  return result;
}

template< class g_func_t, class... Args >
double MassTreeBWPD( MassTree const& mass_tree, g_func_t g_func, Args... args )
{
  /*  The goal is to iterate over all edges and calculate D_s(i), which is
        the fraction of reads in sample s that are in leaves (or edges / placements
        in our case) on the distal side of edge i. Then, the BWPD is the sum of edge
        length l(i) multiplied by [2min(D(i),1−D(i))]^θ. For θ=0, and evaluating only
        those edges that are in a samples spanning tree, this results in the somewhat
        classical FaithPD, except adapted to phylogenetic placement.
    */

  // assert( theta >= 0.0 and theta <= 1.0 );

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
    if( it.is_last_iteration() ) {
      continue;
    }
    ++count;

    auto const& edge       = it.edge();
    auto const& edge_index = edge.index();

    // if this is a leaf, there cannot be any mass on the distal side, so we skip
    if( is_leaf( edge ) ) {
      per_edge_D[ edge_index ] = 0.0;
      continue;
    }

    // interior edges:
    // calculate the new distal count as the sum of distal counts and query counts
    // of the child edges
    auto const& node          = it.node();
    auto const& lhs_edge      = node.link().next().edge();
    auto const& rhs_edge      = node.link().next().next().edge();
    auto const lhs_edge_index = lhs_edge.index();
    auto const rhs_edge_index = rhs_edge.index();

    per_edge_D[ edge_index ] = per_edge_D[ lhs_edge_index ] + mass_per_edge[ lhs_edge_index ] + per_edge_D[ rhs_edge_index ] + mass_per_edge[ rhs_edge_index ];

    auto const branch_length = edge.data< MassTreeEdgeData >().branch_length;

    // update the BWPD sum with the result of this edge
    result += branch_length * g_func( per_edge_D[ edge_index ], args... );
  }

  if( count != mass_tree.edge_count() ) {
    throw std::runtime_error( "count != mass_tree.edge_count()" );
  }

  return result;
}

/**
 *  Calculate different diversity metrics based on the given jplace files, output as CSV
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-files...>\n" );
  }

  std::vector< std::string > jplace_files;
  for( int i = 1; i < argc; ++i ) {
    jplace_files.emplace_back( argv[ i ] );
  }
  SampleSet samples = JplaceReader().read( from_files( jplace_files ) );

  std::vector< double > theta_set = { 0.0, 0.25, 0.5, 0.75, 1.0 };

  // write the header
  std::cout << "sample,phylo_entropy,quadratic";
  for( auto const theta : theta_set ) {
    std::cout << ",bwpd_" << theta;
  }
  std::cout << "\n";

  for( size_t i = 0; i < samples.size(); ++i ) {
    std::cout << samples.name_at( i );
    auto const& sample = samples[ i ];
    // double const n = num_queries( sample );

    // Phylogenetic Entropy
    std::cout << "," << -BWPD( sample, phylo_entropy_g );

    // Phylogenetic Quadratic Entropy
    std::cout << "," << BWPD( sample, phylo_quad_entropy_g );

    // BWPD
    for( auto const theta : theta_set ) {
      std::cout << "," << BWPD( sample, step_function_g, theta );
    }
    std::cout << "\n";
  }

  // trying with the mass_tree
  // write the header
  std::cout << "sample,phylo_entropy,quadratic";
  for( auto const theta : theta_set ) {
    std::cout << ",bwpd_" << theta;
  }
  std::cout << "\n";

  for( size_t i = 0; i < samples.size(); ++i ) {
    std::cout << samples.name_at( i );
    auto& sample = samples[ i ];
    // double const n = num_queries( sample );

    normalize_weight_ratios( sample );

    auto mass_tree = convert_sample_to_mass_tree( sample, true ).first;
    // MassTreeBWPD( mass_tree, step_function_g, 1.0 );
    // break;

    // Phylogenetic Entropy
    std::cout << "," << -MassTreeBWPD( mass_tree, phylo_entropy_g );

    // Phylogenetic Quadratic Entropy
    std::cout << "," << MassTreeBWPD( mass_tree, phylo_quad_entropy_g );

    // BWPD
    for( auto const theta : theta_set ) {
      std::cout << "," << MassTreeBWPD( mass_tree, step_function_g, theta );
    }
    std::cout << "\n";
  }

  return 0;
}
