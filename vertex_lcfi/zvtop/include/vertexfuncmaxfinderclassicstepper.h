#ifndef VERTEXFUNCMAXFINDERDUMBSTEPPER_H
#define VERTEXFUNCMAXFINDERDUMBSTEPPER_H

#include "vertexfuncmaxfinder.h"
#include "../../util/inc/vector3.h"

namespace vertex_lcfi
{
namespace ZVTOP
{
	class VertexFunction;

//!Robust VertexFuncMaxFinder
/*!
Similar to original SLD implementaion, minimises along each axis in turn
by stepping till minimum reached. Currently hard wired step size.
 \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
 \version 0.1
 \date    20/09/05
*/
	class VertexFuncMaxFinderClassicStepper :
		public VertexFuncMaxFinder
	{
	public:
		VertexFuncMaxFinderClassicStepper();
		Vector3 findNearestMaximum(const Vector3 & StartPoint, VertexFunction* VertexFunction);
		VertexFuncMaxFinderClassicStepper(const vertex_lcfi::ZVTOP::VertexFuncMaxFinderClassicStepper&) = delete;
		VertexFuncMaxFinderClassicStepper& operator=(const vertex_lcfi::ZVTOP::VertexFuncMaxFinderClassicStepper&) = delete;
	private:
		
		Vector3 _CurrentPos{};
		double _CurrentValue=0;
		VertexFunction* _Function=nullptr;


		void _minimiseAlongAxis(const Vector3 & Step);
	};
}
}
#endif //VERTEXFUNCMAXFINDERDUMBSTEPPER_H

