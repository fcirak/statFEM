//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2012

//------------------------------------------------------------------------------
/** Insert coordinates into an insert iterator
 *  \tparam IT               Type of the iterator used here
 *  \param[in]      backIter The iterator of some coordinate storage
 *  \return                  The manipulated iterator
 */
//template<typename IT>
//IT gshell::fem::NodeShell::passCoordinates( IT backIter ) const
//{
//    return std::copy( coord_.begin(), coord_.end(), backIter );
//}

//------------------------------------------------------------------------------
/** Write only the ID of this node to a given output stream
 *  \param[in,out]  out Output stream
 *  \return             Output stream
*/
std::ostream & gshell::fem::NodeShell::writeId( std::ostream & out ) const
{
    out << id_ << " ";
    return out;
}

//------------------------------------------------------------------------------
/** Write the Coordinates of this node to a given output stream
 *  \param[in,out]  out Output stream
 *  \return             Output stream
 */
std::ostream & gshell::fem::NodeShell::writeCoordinates( std::ostream & out ) const
{
//    std::ostream_iterator<double> doubleOut( out, "  " );
//    std::copy( coord_.begin(), coord_.end(), doubleOut );
    out << coord_;
    if ( dim == 2 ) out << "  0"; // add another zero for format reasons
    if ( dim == 1 ) out << "  0  0"; // dito
    out << std::endl;
    return out;
}

//------------------------------------------------------------------------------
/** Read in the coordinates from a given input stream
 *  \param[in,out]  inp Input stream
 *  \return             Input stream
 */
std::istream & gshell::fem::NodeShell::readSelf( std::istream & inp )
{
    // read dim coordinates
    for ( unsigned i = 0; i < dim; i++ ) { 
        inp >> coord_[ i ];
        if( inp.peek( ) == ',' ) inp.ignore( );
    }
    // swallow 3-dim coordinates
    for ( unsigned i = dim; i < 3; i++ ) {
        double dummy; 
        inp >> dummy;  
    } 
    return inp;
}

//------------------------------------------------------------------------------
//! Set coefficients of tangent
void gshell::fem::NodeShell::addTangentCoefficient( NodeShell * node,
                                                    const unsigned dir,
                                                    const double coeff )
{
    tangCoeffs_[ dir ][ node ] = coeff;
    return;
}

//------------------------------------------------------------------------------
double gshell::fem::NodeShell::getTangentCoefficient( NodeShell * node,
                                                      const unsigned dir ) const
{
    const MapNPtrDouble & tangCoeff = tangCoeffs_[ dir ];
    MapNPtrDouble::const_iterator nc = tangCoeff.find( node );
    if ( nc != tangCoeff.end() )
        return nc->second;
    else
        return 0.0;
}

//------------------------------------------------------------------------------
//! Set coefficients of limit surface
void gshell::fem::NodeShell::addLimitCoefficient( NodeShell * node,
                                                  const double coeff )
{
    limitCoeff_[ node ] = coeff;
    return;
}

//------------------------------------------------------------------------------
double gshell::fem::NodeShell::getLimitCoefficient( NodeShell * node ) const
{
    MapNPtrDouble::const_iterator nc = limitCoeff_.find( node );
    if ( nc != limitCoeff_.end() )
        return nc->second;
    else
        return 0.0;
}


