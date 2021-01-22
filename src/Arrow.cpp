// breathalyser
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "gArrow.h"
#include "vArrow.h"
#include "fArrow.h"
#include "Arrow.h"
#include "Segment.h"

Arrow::Arrow(Segment *parent)
{
	_segment = parent;
	_renderType = GL_LINES;
	_gString = Arrow_gsh();
	_fString = Arrow_fsh();
	_vString = Arrow_vsh();
	this->SlipObject::setName("Arrow");
}

void Arrow::populate()
{
	vec3 start, dir;
	_segment->motionComparedTo(_segment->sister(), &start, &dir);

	Helen3D::Vertex v;
	memset(&v, '\0', sizeof(Helen3D::Vertex));
	v.color[3] = 1.;

	pos_from_vec(&v.pos[0], start);
	_vertices.push_back(v);

	pos_from_vec(&v.pos[0], dir);
	_vertices.push_back(v);
	
	_indices.push_back(0);
	_indices.push_back(1);
}


