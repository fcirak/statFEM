// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SmfHead.ipp

//------------------------------------------------------------------------------
//! Standrad constructor
corlib::SmfHead::SmfHead( ) : 
    headChar_( '!' ),
    commentChar_( '#' ),
    foundHeader_( false ),
    elementShapeName_( "elementShape" ),
    elementShape_( convertShapeEnumToString( UNDEFINED ) ),
    elementNumPointsName_( "elementNumPoints" ),
    elementNumPoints_( 0 )
{
    parser_.registerPropertiesVar( elementShapeName_, elementShape_ );
    parser_.registerPropertiesVar( elementNumPointsName_, elementNumPoints_ );
    return;
}

//------------------------------------------------------------------------------
//! Try to read the header
std::istream & corlib::SmfHead::read( std::istream & is )
{
    // extract head, 
    // i.e. create an istream in which all '!' are removed, 
    // which appear as first character on a line
    std::stringstream head;
    // a direct cast would be possible as well
    std::istream::int_type i = is.peek();
    FTL_VERIFY( i != is.eof() ); // check state
    char c = static_cast<char>( i ); 

    while ( ( c == headChar_ ) or ( c == commentChar_ ) ) {
        std::string line;
        std::getline( is, line );
        if ( line.at( 0 ) == headChar_ ) {
            foundHeader_ = true;
            line.erase( 0, 1 );
        }
        head << line << '\n';

        std::istream::int_type i = is.peek();
        FTL_VERIFY( i != is.eof() );
        c = static_cast<char>( i ); 
    }
    
    // extract variables
    if ( foundHeader_ ) {
        
        // try to read registered variables (skipping comments)
        parser_.readValues( head );
        
        if ( parser_.hasUnrecognized() ) {
            std::cerr << "Found unrecognized variable(s): " << std::endl;
            parser_.printUnrecognized( std::cerr );
            exit(0);
        }
        
        // check variables have been found
        this->foundVariables();    
    }
    
    // return input stream with skimmed of header
    return is;
}

//------------------------------------------------------------------------------
//! Read and validate quantities
std::istream & corlib::SmfHead::readValidated( std::istream & is,
                                               const enum shape se,
                                               const unsigned np )
{
    // read the header
    this->read( is );
    
    // validate
    if ( foundHeader_ ) {
        FTL_VERIFY_DESCRIPTIVE( ( this->giveElementShape() == se ),
                                "Shape in input file does not match to compiled element shape, i.e. %s != %s\n",
                                corlib::convertShapeEnumToString( this->giveElementShape() ).c_str( ),
                                corlib::convertShapeEnumToString( se ).c_str( ) );
        FTL_VERIFY_DESCRIPTIVE( ( elementNumPoints_ == np ),
                                "Number of element points in input file does not match to compiled number, i.e. %d != %d\n",
                                elementNumPoints_, np );
    }
    
    // return
    return is;
}

//------------------------------------------------------------------------------
//! Check if registered variables have been found
bool corlib::SmfHead::foundVariables( ) const
{
    // initialise
    bool foundVar = true;

    // shape
    if ( elementShape_ == convertShapeEnumToString( UNDEFINED ) ) {
        foundVar = false;
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Element shape was not found in header, though header provided" );
    }
    else {
        foundVar = foundVar and true;
    }

    // number of element points
    if ( elementNumPoints_ == 0 ) {
        foundVar = false;
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Number of element points was not found in header, though header provided" );
    }
    else {
        foundVar = foundVar and true;
    }

    // provide result
    return foundVar;
}

//------------------------------------------------------------------------------
//! write header to stream
void corlib::SmfHead::write( std::ostream & os ) const
{
    this->write( elementShape_, elementNumPoints_, os );
    return;
}

//------------------------------------------------------------------------------
//! write header to stream
void corlib::SmfHead::write( const std::string se,
                             const unsigned np,
                             std::ostream & os ) const
{
    // write header
    os << headChar_ << " " << elementShapeName_ << " " << se << std::endl;
    os << headChar_ << " " << elementNumPointsName_ << " " << np << std::endl;
    
    // done
    return;
}

//------------------------------------------------------------------------------
//! write header to stream
void corlib::SmfHead::write( const enum shape se,
                             const unsigned np,
                             std::ostream & os ) const
{
    this->write( convertShapeEnumToString( se ), np, os );
    return;
}
