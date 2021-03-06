//#include <sdtlib.h>
#include <iostream> //library for C++
#include <QList>	/**<Qt containers*/
#include <QThread>  /**<Qt multicore*/

#include <omp.h>	// OpenMP
#define NTHREADS 2	// Number of CPUs for one job
#include <time.h>

#include <cstdlib>

#include <sys/resource.h>

using namespace std; // for a easier use of cin and cout

#include <math.h>     /**<used for pow(),sqrt()*/
#include <sys/stat.h> /**<used for mkdir();*/

#include <dirent.h>   /**<used for opendir();*/
#include <string> /**<Used to handle string.*/
#include <fstream> /**<Used to open files.*/ 
#include <sstream> /**<Used to convert integer into string.*/

#define DIM 3  /**< Space dimension. */ 
#define PI 3.1415926535 /** Guess what*/
#define MESH_SIZE 1.05  /** ratio of mesh size and the maximum particle diameter */
          
#define PSEUDO_2D 0 /**< 1= yes 0 = no (then 3D)*/	

//AK Addition
#define MAX_CAP_LENGTH 0.5 /**<maximum capillary bridge length before reupture*/
bool MELT_SURFTEN /**<include effect of melt surface tension*/
#define CONTACT_ANGLE_MIN 5.0/180.0*PI /**<minimum contact angle for surface tension*/ 
#define CONTACT_ANGLE_MAX 60.0/180.0*PI /**<maximum contact angle for surface tension*/
#define MELT_SURFACE_TENSION 0.001 /**<surface tension of melt material*/        
	     
class Cvector;
class Cmatrix;
class Cwall;
class Cstat;
class Cpair_mobility_line;
class Cparticle;
class Ccontact;
class Cconfig;
class Cin_out;
// class Caggregate; //YG

class Cparameter;
class Cevent;
class Cpost_process;
class Crun;

class Ccell;
class Cgrid;
class Cnode;
class Cbox;
class Cmesh;
class Cprofile;
#include "matrix.h"
#include "inout.h"
#include "stat.h"
#include "contact.h"
#include "particle.h"
#include "parameter.h"
#include "cell.h"   
#include "mesh.h" 
#include "config.h"
#include "profile.h"
// #include "aggregate.h" //YG

//#include "post_process/grid.h"
//#include "post_process/node.h"
//#include "post_process/distribution.h"
//#include "post_process/post_process.h"	     
#include "run.h"
#include "post_process.h"

#include "matrix.cpp"
#include "stat.cpp"
#include "contact.cpp"
#include "particle.cpp"
#include "inout.cpp"
#include "parameter.cpp"  
#include "cell.cpp"
// #include "aggregate.cpp" //YG
#include "mesh.cpp"
#include "config.cpp"
#include "profile.cpp"
#include "run.cpp"
#include "post_process.cpp"
//#include "post_process/spatial.cpp"
//#include "post_process/grid.cpp"
//#include "post_process/node.cpp"
//#include "post_process/distribution.cpp"
//#include "post_process/post_process.cpp"	    
	
	
	/*! \mainpage 

\author Pierre Rognon
\version thermal-dry
\date    May-2010

 \section intro_sec Introduction
 This version of the Soft-Dynamics can simulate the motion and the temperature of dry grains, including:

  \li  elastic-frictional contact, normal viscous dissipation, rolling and twist torque and friction;
  \li  heat conduction through contacts;
  \li  shear flow;
  \li  flow down a slope under gravity.

\section licence_sec Licence

The Soft_Dynamics is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
   
    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License in the file 'COPYING' along with the Soft Dynamics package. If not, see <http://www.gnu.org/licenses/>.

\section install_sec Installation-compilation

The Soft Dynamics is written in C++ and include other GPL library (see file soft_dynamics.h). It uses the Qt library, and can be easily complile using the qmake facility, for instance through Qdevelop editor.
*/     
