/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "enumerate_linear_transformations.hpp"

#include <fstream>
#include <iostream>
#include <unordered_map>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <igraph/igraph.h>

#include <core/utils/bitset_utils.hpp>
#include <core/utils/range_utils.hpp>

#define L(x) { if ( verbose ) { std::cout << x << std::endl; } }

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

struct boolean_square_matrix
{
public:
  boolean_square_matrix( const boost::dynamic_bitset<>& matrix, unsigned dimension )
    : matrix( matrix ),
      dimension( dimension )
  {
  }

  boolean_square_matrix( const std::vector<boost::dynamic_bitset<>>& rows )
  {
    assert( !rows.empty() );

    dimension = rows.front().size();
    matrix.resize( dimension * dimension );

    unsigned index = 0u;

    for ( const auto& r : rows )
    {
      for ( auto c = 0u; c < dimension; ++c )
      {
        matrix[index++] = r[c];
      }
    }
  }

  static boolean_square_matrix identity( unsigned dimension )
  {
    std::vector<boost::dynamic_bitset<>> rows( dimension, boost::dynamic_bitset<>( dimension ) );
    for ( auto c = 0u; c < dimension; ++c )
    {
      rows[c].set( c );
    }
    return boolean_square_matrix( rows );
  }

  inline bool at( unsigned row, unsigned column ) const
  {
    return matrix.test( row * dimension + column );
  }

  inline void set( unsigned row, unsigned column, bool value )
  {
    matrix.set( row * dimension + column, value );
  }

  std::vector<boost::dynamic_bitset<>> get_row_vectors() const
  {
    std::vector<boost::dynamic_bitset<>> rows( dimension, boost::dynamic_bitset<>( dimension ) );

    for ( auto r = 0u; r < dimension; ++r )
    {
      for ( auto c = 0u; c < dimension; ++c )
      {
        rows[r][c] = at( r, c );
      }
    }

    return rows;
  }

  boolean_square_matrix get_submatrix_by_removing( unsigned row, unsigned column ) const
  {
    boost::dynamic_bitset<> subm( ( dimension - 1u ) * ( dimension - 1u ) );

    unsigned index = 0;
    for ( auto r = 0u; r < dimension; ++r )
    {
      if ( r == row ) continue;
      for ( auto c = 0u; c < dimension; ++c )
      {
        if ( c == column ) continue;
        subm[index++] = at( r, c );
      }
    }

    return boolean_square_matrix( subm, dimension - 1u );
  }

  bool determinant()
  {
    if ( dimension == 2u )
    {
      return ( at( 0, 0 ) && at( 1, 1 ) ) != ( at( 0, 1 ) && at( 1, 0 ) );
    }
    else
    {
      bool result = false;
      for ( auto c = 0u; c < dimension; ++c )
      {
        result = result != ( at( 0u, c ) * get_submatrix_by_removing( 0u, c ).determinant() );
      }
      return result;
    }
  }

  bool is_sorted()
  {
    const auto rows = get_row_vectors();

    auto v = 0u;

    for ( const auto& r : rows )
    {
      const auto rv = r.to_ulong();
      if ( rv < v ) return false;
      v = rv;
    }

    return true;
  }

  boolean_square_matrix swap( unsigned row1, unsigned row2 ) const
  {
    boolean_square_matrix other = *this;

    for ( auto c = 0u; c < dimension; ++c )
    {
      other.set( row1, c, at( row2, c ) );
      other.set( row2, c, at( row1, c ) );
    }

    return other;
  }

  boolean_square_matrix add_to( unsigned from_row, unsigned to_row ) const
  {
    boolean_square_matrix other = *this;

    for ( auto c = 0u; c < dimension; ++c )
    {
      other.set( to_row, c, at( from_row, c ) != at( to_row, c ) );
    }

    return other;
  }

  boolean_square_matrix sort() const
  {
    auto rows = get_row_vectors();

    std::sort( rows.begin(), rows.end(), []( const boost::dynamic_bitset<>& bs1, const boost::dynamic_bitset<>& bs2 ) { return bs1.to_ulong() < bs2.to_ulong(); } );

    return boolean_square_matrix( rows );
  }

  /* multiply a column vector */
  boost::dynamic_bitset<> multiply( const boost::dynamic_bitset<>& v ) const
  {
    boost::dynamic_bitset<> y( dimension );

    for ( auto r = 0u; r < dimension; ++r )
    {
      auto value = false;
      for ( auto c = 0u; c < dimension; ++c )
      {
        value = value != ( v[c] && at( r, c ) );
      }
      y[r] = value;
    }

    return y;
  }                                

  std::pair<unsigned, unsigned> compute_swap_permutation( const boolean_square_matrix& other ) const
  {
    std::vector<unsigned> perm( 1u << dimension ); /* permutation of this matrix */
    std::vector<unsigned> perm2( 1u << dimension );

    boost::dynamic_bitset<> bs( dimension );

    do
    {
      perm[bs.to_ulong()] = multiply( bs ).to_ulong();
      
      inc( bs );
    } while ( bs.any() );

    do
    {
      perm2[other.multiply( bs ).to_ulong()] = bs.to_ulong();

      inc( bs );
    } while ( bs.any() );

    unsigned delta = 0u;
    unsigned theta = 0u;

    for ( auto i = 0u; i < perm.size(); ++i )
    {
      const auto j = perm2[perm[i]];
      if ( i < j )
      {
        delta = j - i;
        theta |= ( 1 << i );
      }
    }

    return {delta, theta};
  }

  inline unsigned long value() const
  {
    return matrix.to_ulong();
  }

private:
  boost::dynamic_bitset<> matrix;
  unsigned                dimension;
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

unsigned factorial( unsigned n )
{
  return ( n == 1u || n == 0u ) ? 1u : factorial( n - 1u ) * n;
}

unsigned compute_order( unsigned n )
{
  unsigned order = 1u;

  for ( auto i = 0; i < n; ++i )
  {
    order *= ( 1 << n ) - ( 1 << i );
  }

  return order;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void enumerate_linear_transformations( unsigned n, const properties::ptr& settings, const properties::ptr& statistics )
{
  /* settings */
  const auto abstract = get( settings, "abstract", false );
  const auto dotname  = get( settings, "dotname",  std::string() );
  const auto verbose  = get( settings, "verbose",  false );
  
  std::vector<boolean_square_matrix> node_to_matrix;
  std::unordered_map<unsigned long, unsigned> matrix_to_node;

  igraph_t graph, * p_g = &graph;

  auto order = compute_order( n );
  if ( abstract )
  {
    order /= factorial( n );
  }
  L( "[i] order = " << order );
  igraph_empty( p_g, order, IGRAPH_UNDIRECTED );
  
  boost::dynamic_bitset<> bs( n * n );

  /* create vertices */
  foreach_bitset( n * n, [&matrix_to_node, &node_to_matrix, n, abstract, verbose]( const boost::dynamic_bitset<>& bs ) {
      boolean_square_matrix m( bs, n );
      if ( m.determinant() != 0u && ( !abstract || m.is_sorted() ) )
      {
        matrix_to_node[m.value()] = node_to_matrix.size();
        node_to_matrix.push_back( m );
      }
    } );

  /* create edges */
  for ( auto v = 0; v < node_to_matrix.size(); ++v )
  {
    const auto& bm = node_to_matrix[v];

    /* adding */
    for ( auto from_row = 0u; from_row < n; ++from_row )
    {
      for ( auto to_row = 0u; to_row < n; ++to_row )
      {
        if ( from_row == to_row ) continue;

        auto other = bm.add_to( from_row, to_row );
        if ( abstract )
        {
          other = other.sort();
        }

        const auto w = matrix_to_node[other.value()];
        if ( w > v )
        {
          igraph_add_edge( p_g, v, w );
        }
      }
    }

    /* swapping */
    if ( abstract )
    {
      for ( auto row2 = 1; row2 < n; ++row2 )
      {
        for ( auto row1 = 0; row1 < row2; ++row1 )
        {
          const auto other = bm.swap( row1, row2 );
          const auto w = matrix_to_node[other.value()];
          if ( w > v )
          {
            igraph_add_edge( p_g, v, w );
          }
        }
      }
    }
  }

  /* check some properties */
  igraph_bool_t b_connected;
  igraph_is_connected( p_g, &b_connected, IGRAPH_STRONG );
  L( "[i] is strongly connected: " << b_connected );

  /* Hamiltonian path */
  igraph_t rgraph, * p_rg = &rgraph;
  igraph_ring( p_rg, order, 0, 0, 1 );

  igraph_bool_t b_iso;
  igraph_vector_t map;
  igraph_vector_init( &map, 0 );
  igraph_subisomorphic_lad( p_rg, p_g, nullptr, &b_iso, &map, nullptr, 0, 0 );
  L( "[i] found isomorphism: " << b_iso );

  igraph_vector_print( &map );
  igraph_vector_destroy( &map );
  igraph_destroy( p_rg );

  if ( !dotname.empty() )
  {
    auto * fp = fopen( dotname.c_str(), "w" );
    igraph_write_graph_dot( p_g, fp );
    fclose( fp );
  }

  igraph_destroy( p_g );
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
