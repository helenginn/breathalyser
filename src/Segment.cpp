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

#include <libsrc/Atom.h>
#include "Segment.h"
#include <h3dsrc/Mesh.h>
#include <QThread>

Segment::Segment() : QObject(), SlipObject(), AtomGroup()
{
	_w = NULL;
}

void Segment::startEnd(int *min, int *max)
{
	*min = INT_MAX;
	*max = -INT_MAX;

	for (size_t i = 0; i < atomCount(); i++)
	{
		int resi = atom(i)->getResidueNum();
		if (resi > *max)
		{
			*max = resi;
		}
		
		if (resi < *min)
		{
			*min = resi;
		}
	}
}

void Segment::populate()
{
	vec3 c = AtomGroup::centroid();

	for (size_t i = 0; i < atomCount(); i++)
	{
		vec3 a = atom(i)->getAbsolutePosition();
		vec3 diff = vec3_subtract_vec3(a, c);
		vec3_set_length(&diff, 1);
		
		Helen3D::Vertex v;
		memset(&v, '\0', sizeof(Helen3D::Vertex));
		pos_from_vec(&v.pos[0], a);
		pos_from_vec(&v.normal[0], diff);
		
		_vertices.push_back(v);
	}
	
	/* you must now make a mesh and link populate() to other things */
	
	makeMesh();
	refineMesh();
}

void Segment::refineMesh()
{
	if (_w == NULL)
	{
		_w = new QThread();
		mesh()->moveToThread(_w);
	}
	
	_w->start();
	connect(mesh(), SIGNAL(resultReady()), this, SLOT(handleMesh()));
	connect(this, SIGNAL(refine()), mesh(), SLOT(shrinkWrap()));
	
	emit refine();
}

void Segment::handleMesh()
{
	disconnect(this, SIGNAL(refine()), nullptr, nullptr);
	disconnect(nullptr, nullptr, this, SLOT(handleMesh()));
	_w->quit();
	
	mesh()->changeToTriangles();
}

void Segment::render(SlipGL *gl)
{
	mesh()->render(gl);
}

