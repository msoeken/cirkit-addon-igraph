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
#include <boost/format.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <igraph/igraph.h>

#include <core/utils/bitset_utils.hpp>
#include <core/utils/graph_utils.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/string_utils.hpp>
#include <classical/utils/truth_table_utils.hpp>

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

  unsigned row_sum() const
  {
    boost::dynamic_bitset<> sum( dimension );
    for ( const auto& r : get_row_vectors() )
    {
      for ( auto c = 0u; c < dimension; ++c )
      {
        sum[c] = sum[c] != r[c];
      }
    }
    return sum.to_ulong();
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

  /* other can be sorted, but this must be as we want it
     we assume it's possible */
  std::pair<unsigned, unsigned> get_add_to_rows( const boolean_square_matrix& other ) const
  {
    const auto vectors = get_row_vectors();
    const auto ovectors = other.get_row_vectors();

    auto to = 0u, from = 0u;
    boost::dynamic_bitset<> to_row;
    boost::dynamic_bitset<> from_row;

    for ( const auto& r : index( vectors ) )
    {
      const auto it = std::find( ovectors.begin(), ovectors.end(), r.value );
      if ( it == ovectors.end() )
      {
        to = r.index;
        to_row = r.value;
      }
    }

    for ( const auto& r : ovectors )
    {
      const auto it = std::find( vectors.begin(), vectors.end(), r );
      if ( it == vectors.end() )
      {
        from_row = r;
      }
    }

    from_row ^= to_row;
    from = std::distance( vectors.begin(), std::find( vectors.begin(), vectors.end(), from_row ) );

    //std::cout << boost::format( "[i] added %s [%d] to %s [%d]" ) % to_string( from_row ) % from % to_string( to_row ) % to << std::endl;

    return {from, to};
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

  friend std::ostream& operator<<( std::ostream& os, const boolean_square_matrix& m )
  {
    for ( auto r = 0u; r < m.dimension; ++r )
    {
      for ( auto c = 0u; c < m.dimension; ++c )
      {
        if ( c > 0u )
        {
          os << " ";
        }
        os << m.at( r, m.dimension - c - 1 );
      }
      os << std::endl;
    }

    return os;
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

  for ( unsigned i = 0; i < n; ++i )
  {
    order *= ( 1 << n ) - ( 1 << i );
  }

  return order;
}

int igraph_is_strongly_regular( const igraph_t *graph, igraph_bool_t *res, igraph_integer_t *k, igraph_integer_t *lambda, igraph_integer_t *mu )
{
  int err;
  *lambda = -1;
  *mu = -1;

  /* check g is undirected */
  if ( igraph_is_directed( graph ) ) { return IGRAPH_EINVAL; }

  /* special case if empty */
  if ( igraph_vcount( graph ) == 0u )
  {
    *res = 1;
    return IGRAPH_SUCCESS;
  }

  /* check same degree property */
  igraph_vector_t degrees;
  igraph_vector_init( &degrees, 0 );
  igraph_vs_t vertices;
  igraph_vs_all( &vertices );
  if ( ( err = igraph_degree( graph, &degrees, vertices, IGRAPH_ALL, 0 ) ) != IGRAPH_SUCCESS )
  {
    //return err;
  }
  igraph_vs_destroy( &vertices );

  *k = VECTOR( degrees )[0];
  for ( auto i = 1u; i < igraph_vector_size( &degrees ); ++i )
  {
    if ( VECTOR( degrees )[i] != *k )
    {
      igraph_vector_destroy( &degrees );
      std::cout << "[w] not all vertices have same degree" << std::endl;
      *res = 0;
      return IGRAPH_SUCCESS;
    }
  }

  igraph_vector_destroy( &degrees );

  /* check that neighbors have the same number of neighbors */
  for ( auto node = 0u; node < igraph_vcount( graph ); ++node )
  {
    igraph_vector_t neis;
    igraph_vector_init( &neis, 0 );
    if ( ( err = igraph_neighbors( graph, &neis, node, IGRAPH_ALL ) ) != IGRAPH_SUCCESS )
    {
      // return err;
    }

    /* point to the first neighbor */
    auto cur = 0;
    auto conn = 0;
    while ( cur < igraph_vector_size( &neis ) && VECTOR( neis )[cur] < node + 1 ) ++cur;

    for ( auto onode = node + 1; onode < igraph_vcount( graph ); ++onode )
    {
      if ( onode == node ) continue;

      if ( cur < igraph_vector_size( &neis ) && VECTOR( neis )[cur] == onode )
      {
        /* connected */
        conn = 1;
        ++cur;
      }
      else
      {
        /* not connected */
        conn = 0;
      }

      igraph_vector_t oneis;
      igraph_vector_init( &oneis, 0 );
      igraph_neighbors( graph, &oneis, onode, IGRAPH_ALL );

      igraph_vector_t inter;
      igraph_vector_init( &inter, 0 );
      igraph_vector_intersect_sorted( &neis, &oneis, &inter );

      if ( conn )
      {
        if ( *lambda == -1 )
        {
          *lambda = igraph_vector_size( &inter );
        }
        else if ( *lambda != igraph_vector_size( &inter ) )
        {
          igraph_vector_destroy( &inter );
          igraph_vector_destroy( &oneis );
          igraph_vector_destroy( &neis );

          std::cout << "[w] not same lambda" << std::endl;
          *res = 0;
          return IGRAPH_SUCCESS;
        }
      }
      else
      {
        if ( *mu == -1 )
        {
          *mu = igraph_vector_size( &inter );
        }
        else if ( *mu != igraph_vector_size( &inter ) )
        {
          igraph_vector_destroy( &inter );
          igraph_vector_destroy( &oneis );
          igraph_vector_destroy( &neis );

          std::cout << "[w] not same mu" << std::endl;
          *res = 0;
          return IGRAPH_SUCCESS;
        }
      }

      igraph_vector_destroy( &inter );

      igraph_vector_destroy( &oneis );
    }

    igraph_vector_destroy( &neis );
  }

  *res = 1;
  return IGRAPH_SUCCESS;
}

struct enumerate_hamiltonians_is_equivalent
{
  enumerate_hamiltonians_is_equivalent( const std::vector<boolean_square_matrix>& node_to_matrix )
    : node_to_matrix( node_to_matrix )
  {
  }

  bool operator()( const vertex_t<graph_t<>>& vs, const vertex_t<graph_t<>>& vl ) const
  {
    if ( vs == 0u )
    {
      return vl == 0u;
    }
    if ( vs < prefix.size() && prefix[vs] != 0 )
    {
      return node_to_matrix[vl].row_sum() == prefix[vs];
    }
    else
    {
      return true;
    }
  }

private:
  //const std::vector<unsigned> prefix = {1, 2, 4, 3, 6, 7, 5, 3, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 7};
  //const std::vector<unsigned> prefix = {1, 2, 4, 3, 6, 7, 5, 3, 2, 4, 7, 1, 5, 6, 5, 4, 6, 1, 3, 7, 2, 7, 4, 1, 3, 5, 2, 6};
  const std::vector<unsigned> prefix = {7, 6, 5, 2, 4, 1, 3, 5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3};
  const std::vector<boolean_square_matrix>& node_to_matrix;
};

struct enumerate_hamiltonians_printer
{
public:
  enumerate_hamiltonians_printer( const graph_t<>& sg, const std::vector<boolean_square_matrix>& node_to_matrix )
    : sg( sg ),
      node_to_matrix( node_to_matrix )
  {
  }

  template<typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()( const CorrespondenceMap1To2& f, const CorrespondenceMap2To1& g ) const
  {
    std::vector<std::vector<unsigned>> parts( 4u, std::vector<unsigned>( 7u ) );

    for ( auto v = 0u; v < boost::num_vertices( sg ); ++v )
    {
      parts[v / 7][v % 7] = node_to_matrix[boost::get( f, v )].row_sum();
    }

    for ( auto& p : parts )
    {
      std::sort( p.begin(), p.end() );
      for ( auto j = 0u; j < 7u; ++j )
      {
        if ( p[j] != j + 1 )
        {
          return true;
        }
      }
    }

    for ( auto v = 0u; v < boost::num_vertices( sg ); ++v )
    {
      std::cout << " " << node_to_matrix[boost::get( f, v )].row_sum();
    }
    std::cout << std::endl;

    for ( auto v = 0u; v < boost::num_vertices( sg ); ++v )
    {
      std::cout << " " << boost::get( f, v );
    }
    std::cout << std::endl;

    for ( auto i = 0u; i < 4u; ++i )
    {
      std::cout << " " << node_to_matrix[boost::get( f, 7 * i )].row_sum();
    }
    std::cout << std::endl << std::endl;

    return true;
  }

private:
  const graph_t<>& sg;
  const std::vector<boolean_square_matrix>& node_to_matrix;
};

void enumerate_hamiltonians( const igraph_t* graph, const std::vector<boolean_square_matrix>& node_to_matrix )
{
  graph_t<> g;
  for ( auto i = 0; i < igraph_vcount( graph ); ++i )
  {
    add_vertex( g );
  }

  igraph_es_t es;
  igraph_es_all( &es, IGRAPH_EDGEORDER_ID );

  igraph_eit_t eit;
  igraph_eit_create( graph, es, &eit );

  while ( !IGRAPH_EIT_END( eit ) )
  {
    const auto eid = IGRAPH_EIT_GET( eit );

    igraph_integer_t from, to;
    igraph_edge( graph, eid, &from, &to );

    add_edge( from, to, g );

    IGRAPH_EIT_NEXT( eit );
  }

  igraph_eit_destroy( &eit );
  igraph_es_destroy( &es );

  const auto rg = ring_graph( boost::num_vertices( g ) );
  boost::vf2_subgraph_mono( rg, g,
                            enumerate_hamiltonians_printer( rg, node_to_matrix ),
                            boost::vertex_order_by_mult( rg ),
                            boost::vertices_equivalent( enumerate_hamiltonians_is_equivalent( node_to_matrix ) ) );
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void enumerate_linear_transformations( unsigned n, const properties::ptr& settings, const properties::ptr& statistics )
{
  /* settings */
  const auto abstract    = get( settings, "abstract",    true );
  const auto hamiltonian = get( settings, "hamiltonian", false );
  const auto hamilenum   = get( settings, "hamilenum",   false );
  const auto frompath    = get( settings, "frompath",    std::string() );
  const auto dotname     = get( settings, "dotname",     std::string() );
  const auto verbose     = get( settings, "verbose",     false );

  std::vector<boolean_square_matrix> node_to_matrix;
  std::unordered_map<unsigned long, unsigned> matrix_to_node;
  std::vector<unsigned> edge_to_operation;

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
  for ( size_t v = 0; v < node_to_matrix.size(); ++v )
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

          igraph_integer_t eid;
          igraph_get_eid( p_g, &eid, v, w, 0, 0 );
          edge_to_operation.push_back( from_row * 3 + to_row );
        }
      }
    }

    /* swapping */
    if ( !abstract )
    {
      for ( unsigned row2 = 1u; row2 < n; ++row2 )
      {
        for ( unsigned row1 = 0u; row1 < row2; ++row1 )
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

  igraph_bool_t b_strongly_regular;
  igraph_integer_t k, lambda, mu;
  assert( !igraph_is_strongly_regular( p_g, &b_strongly_regular, &k, &lambda, &mu ) );
  L( "[i] is strongly regular: " << b_strongly_regular );

  /* Hamiltonian path */
  if ( hamiltonian )
  {
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
  }

  if ( hamilenum )
  {
    enumerate_hamiltonians( p_g, node_to_matrix );
  }

  /* derive delta-swap sequence from path */
  if ( !frompath.empty() )
  {
    std::cout << "[i] derive delta-swap sequence from path: " << frompath << std::endl;
    std::vector<unsigned> path;
    parse_string_list( path, frompath, " " );

    unsigned delta, theta;
    const auto& swaps = tt_store::i().swaps( n );

    auto curm = node_to_matrix[path[0u]];

    std::cout << curm;
    for ( auto i = 0u; i < swaps.size(); ++i )
    {
      auto nextm = curm.swap( swaps[i], swaps[i] + 1 );
      std::tie( delta, theta ) = nextm.compute_swap_permutation( curm );
      std::cout << boost::format( "[i] swap %d with %d (%d, %d)" ) % swaps[i] % ( swaps[i] + 1 ) % delta % theta << std::endl;
      std::cout << nextm;
      curm = nextm;
    }

    auto next = 1u;
    while ( next < path.size() )
    {
      // const auto v = std::min( path[next - 1], path[next] );
      // const auto w = std::max( path[next - 1], path[next] );
      // igraph_integer_t eid;
      // igraph_get_eid( p_g, &eid, v, w, 0, 0 );
      // const auto op = edge_to_operation[eid];
      // const auto to = op % n;
      // const auto from = ( op - to ) / n;
      // std::cout << boost::format( "[i] add %d to %d" ) % from % to << std::endl;

      unsigned from_row, to_row;
      std::tie( from_row, to_row ) = curm.get_add_to_rows( node_to_matrix[path[next]] );

      auto nextm = curm.add_to( from_row, to_row );
      std::tie( delta, theta ) = nextm.compute_swap_permutation( curm );
      std::cout << boost::format( "[i] add %d to %d (%d, %d)" ) % from_row % to_row % delta % theta << std::endl;
      std::cout << nextm;
      curm = nextm;

      for ( auto i = 0u; i < swaps.size(); ++i )
      {
        auto nextm = curm.swap( swaps[i], swaps[i] + 1 );
        std::tie( delta, theta ) = nextm.compute_swap_permutation( curm );
        std::cout << boost::format( "[i] swap %d with %d (%d, %d)" ) % swaps[i] % ( swaps[i] + 1 ) % delta % theta << std::endl;
        std::cout << nextm;
        curm = nextm;
      }
      ++next;
    }
  }

  /* row value inspection */
  if ( false )
  {
    const auto from_sum = 2u;
    const auto to_sum = 3u;
    for ( const auto& bm : node_to_matrix )
    {
      if ( bm.row_sum() == from_sum )
      {
        std::cout << "[i] matrix with row sum " << from_sum << ": " << std::endl << bm;
        //std::cout << "[i] adjacent matrices have row sums:";

        igraph_vector_t neis;
        igraph_vector_init( &neis, 0 );
        igraph_neighbors( p_g, &neis, matrix_to_node[bm.value()], IGRAPH_ALL );

        for ( auto i = 0u; i < igraph_vector_size( &neis ); ++i )
        {
          const auto& om = node_to_matrix[VECTOR( neis )[i]];
          if ( om.row_sum() == to_sum )
          {
            std::cout << "[i] connects to: " << std::endl << om;
          }
        }
        //std::cout << std::endl;
        igraph_vector_destroy( &neis );
      }
    }
  }

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
