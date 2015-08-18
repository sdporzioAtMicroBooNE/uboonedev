/// \file    GeometryDrawer.h
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
/// \version $Id: GeometryDrawer.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
#ifndef EVD_GEOMETRYTPC_H
#define EVD_GEOMETRYTPC_H
#include <vector>

namespace evdb3D
{
    /// Aid in the rendering of Geometry objects
    class GeometryTPC
    {
    public:
        GeometryTPC();
        ~GeometryTPC();
        void DetOutline3D();
    };
}

#endif
////////////////////////////////////////////////////////////////////////
