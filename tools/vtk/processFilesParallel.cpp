#include <iostream>
#include <iterator>
#include <fstream>
#include <mpi.h>

#include <corlib/verify.hpp>

#include <vector>
#include <boost/lexical_cast.hpp>

/* Simple tool to process distributed VTU data files in parallel.
 * Operations are NOT distributed with respect to subdomains but with respect to timesteps. 
 */
int main( int argc, char * argv[] )
{
    if( argc != 6 ) {

        std::cout << "Usage: " << argv[0]  
                  << "  file_with_time_steps" << std::endl
                  << "  basename" << std::endl
                  << "  lower    - lower string of number range " << std::endl
                  << "  upper    - upper string of number range " << std::endl
                  << "  gridfile ( see vtpProbe for details ) " << std::endl;
        return 1;
    }

    // initialize
    MPI_Init( &argc, &argv );
    MPI_Comm commAll = MPI_COMM_WORLD;

    // orient in communicator
    int rank, nProc;
    MPI_Comm_rank( commAll, &rank );
    MPI_Comm_size( commAll, &nProc );


    // load strings with timesteps user wants to convert
    const std::string timeStepsFile = boost::lexical_cast<std::string>( argv[1] );
    const std::string basename      = boost::lexical_cast<std::string>( argv[2] );
    const std::string lower         = boost::lexical_cast<std::string>( argv[3] );
    const std::string upper         = boost::lexical_cast<std::string>( argv[4] );
    const std::string gridfile      = boost::lexical_cast<std::string>( argv[5] );



    std::vector<std::string> timeSteps;

    std::ifstream strts( timeStepsFile.c_str() );
    FTL_VERIFY( strts.is_open() );

    std::string auxStr;
    strts >> auxStr;
    while (!strts.eof( ))     //if not at end of file, continue reading timesteps
    {
        timeSteps.push_back( auxStr );
        strts >> auxStr;
    }
    strts.close( );

    int numSteps = timeSteps.size();
    //std::cout << "numSteps:" << numSteps << std::endl;

    //std::copy( timeSteps.begin(), timeSteps.end(), std::ostream_iterator<std::string>( std::cout, " " ) );
    
    // distribute work
    int numStepsLocMax = ( numSteps + nProc - 1) / nProc;
    int numStepsLoc    = std::min( numSteps - rank * numStepsLocMax, numStepsLocMax );
    int start  = numStepsLocMax * rank;
    int finish = start + numStepsLoc;
    std::cout << "rank:" << rank << "start " << start << "finish" << finish << std::endl;

    std::string command;
    int ierr;
    for ( unsigned iStep = start; iStep < finish; ++iStep ) {

        // merge files for timestep
        command = "${OPENFTLROOT}/tools/vtk/vtuMerge"" " 
                + basename + "."" " 
                + lower + " " 
                + upper + " " 
                + "-" + timeSteps[iStep] + ".vtu";
        std::cout << "Calling system command: " << command << std::endl;
        ierr = system( command.c_str() );
        FTL_VERIFY_DESCRIPTIVE( ierr == 0, "Unable to process merge command." );

        // convert merged VTU file to VTP
        command = "${OPENFTLROOT}/tools/vtk/vtu2vtp"" " 
                + basename + ".merge-" 
                + timeSteps[iStep] + ".vtu";
        std::cout << "Calling system command: " << command << std::endl;
        ierr = system( command.c_str() );
        FTL_VERIFY_DESCRIPTIVE( ierr == 0, "Unable to process vtu2vtp command." );

        // probing required data on rectangular grid
        command = "${OPENFTLROOT}/tools/vtk/vtpProbe"" " 
                + basename + ".merge-" 
                + timeSteps[iStep] + ".vtp"" "
                + gridfile;
        std::cout << "Calling system command: " << command << std::endl;
        ierr = system( command.c_str() );
        FTL_VERIFY_DESCRIPTIVE( ierr == 0, "Unable to process vtpProbe command." );
    }

    MPI_Finalize();

    return 0;
}
