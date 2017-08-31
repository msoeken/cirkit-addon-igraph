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

#include "dswap.hpp"

#include <boost/program_options.hpp>

#include <core/utils/program_options.hpp>

#include <classical/enumeration/enumerate_linear_transformations.hpp>

namespace cirkit
{

using boost::program_options::value;

dswap_command:: dswap_command( const environment::ptr& env )
  : cirkit_command( env, "Find delta-swap cycles" )
{
  opts.add_options()
    ( "num_vars,n", value_with_default( &num_vars ), "number of variables" )
          ("abstract","work on the equivalence classes")
          ("hamiltonian", "check for and compute a hamiltonian path")
          ("hamilenum", "enumerate hamiltonian paths")
    ( "frompath",   value( &frompath ),              "generate delta-swap sequence from path" )
    ( "dotname",    value( &dotname ),               "write graph to this file" )
    ;
  add_positional_option( "num_vars" );
  be_verbose();
}

bool dswap_command::execute()
{
  const auto settings = make_settings();

  if ( is_set( "dotname" ) )
  {
    settings->set( "dotname", dotname );
  }
  if ( is_set( "frompath" ) )
  {
    settings->set( "frompath", frompath );
  }

  settings->set("abstract", is_set("abstract"));
  settings->set("hamiltonian", is_set("hamiltonian"));
  settings->set("hamilenum", is_set("hamilenum"));



  enumerate_linear_transformations( num_vars, settings );

  return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
