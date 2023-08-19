
#ifndef BresenhamIterator_h
#define BresenhamIterator_h

#include <aRibeiroCore/common.h>

namespace aRibeiro {

/// \brief Classical line iterator algorithm
///
/// Generate a raster of a line specified by two points.
///
/// Example:
///
/// \code
/// #include <aRibeiroCore/aRibeiroCore.h>
/// using namespace aRibeiro;
///
/// bool matrix[height][width];
/// int x, y, xd yd;
/// 
/// BresenhamIterator it(x, y, xd, yd);
/// 
/// it.getXY(&x, &y);
/// matrix[y][x] = true;
/// while (it.next(&x, &y)) {
///     matrix[y][x] = true;
/// }
/// \endcode
///
/// \author Alessandro Ribeiro
///
/// https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
class BresenhamIterator{
    
    int dx, dy, sx, sy, error;
    int x0, y0, x1, y1;

    public:
    /// \brief Initialize the Bresenham algorithm state
    ///
    /// It is needed to specify the source and target point as integers.
    ///
    /// There are several algorithms to rasterize lines, but the Bresenham
    /// is eficient and it dont uses floating points to compute the intermediate
    /// points to make the line.
    ///
    /// \author Alessandro Ribeiro
    /// \param x origin x
    /// \param y origin y
    /// \param xd target x
    /// \param yd target y
    ///
    ARIBEIRO_INLINE BresenhamIterator(int _x0, int _y0, int _x1, int _y1) {
        x0 = _x0;
        y0 = _y0;
        x1 = _x1;
        y1 = _y1;
        if (x0 < x1) {
            sx = 1;
            dx = x1 - x0;
        }
        else {
            sx = -1;
            dx = x0 - x1;
        }
        if (y0 < y1) {
            sy = 1;
            dy = y0 - y1;
        }
        else {
            sy = -1;
            dy = y1 - y0;
        }
        error = dx + dy;
    }
    /// \brief Returns true until the algorithm reachs the target point
    /// \author Alessandro Ribeiro
    /// \return True if the target point is reached by the iteration
    ///
    ARIBEIRO_INLINE bool next() {

        if (x0 == x1 && y0 == y1)
            return false;

        int e2 = 2 * error;
        bool r = false;
        if (e2 >= dy) {
            if (x0 == x1) return false;
            error = error + dy;
            x0 = x0 + sx;
            r = true;
        }
        if (e2 <= dx) {
            if (y0 == y1) return false;
            error = error + dx;
            y0 = y0 + sy;
            r = true;
        }

        return r;
    }
    /// \brief Reads the current point of the rasterization state
    /// \author Alessandro Ribeiro
    /// \param x Returns the current x of the iteration
    /// \param y Returns the current y of the iteration
    ///
    ARIBEIRO_INLINE void getXY(int *x,int *y) {
        *x = x0;
        *y = y0;
    }

    /// \brief Returns true until the algorithm reachs the target point.
    /// Reads the current point of the rasterization state
    /// \author Alessandro Ribeiro
    /// \param x Returns the current x of the iteration
    /// \param y Returns the current y of the iteration
    /// \return True if the target point is reached by the iteration
    ///
    ARIBEIRO_INLINE bool next(int *x,int *y) {
        bool retorno = next();
        *x = x0;
        *y = y0;
        return retorno;
    }

};

}

#endif
