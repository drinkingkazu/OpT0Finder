//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//
// Functions
//
#ifndef __CINT__
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< short              > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< unsigned short     > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< int                > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< unsigned int       > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< long long          > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< unsigned long long > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< float              > &)+;
#pragma link C++ function larcv::copy_array(PyObject*, const std::vector< double             > &)+;
//#pragma link C++ function larcv::as_mat(const larcv::Image2D&)+;
//#pragma link C++ function larcv::as_ndarray<short>(const std::vector<short>&)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< short              >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< unsigned short     >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< int                >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< unsigned int       >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< long long          >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< unsigned long long >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< float              >& vec)+;
#pragma link C++ function larcv::as_ndarray (const std::vector< double             >& vec)+;
#endif
//ADD_NEW_CLASS ... do not change this line

#endif




















