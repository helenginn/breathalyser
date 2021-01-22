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
#include <libsrc/Monomer.h>
#include <h3dsrc/Mesh.h>
#include <h3dsrc/SlipGL.h>
#include <hcsrc/maths.h>
#include <QThread>

#include <h3dsrc/shaders/vStructure.h>
#include <h3dsrc/shaders/fStructure.h>

#include "Arrow.h"
#include "Segment.h"
#include "Difference.h"

Segment::Segment(Difference *d, double threshold) 
: QObject(), SlipObject(), AtomGroup()
{
	_w = NULL;
	_threshold = threshold;
	_diff = d;
	_meshDot = 0.5;
	_sister = NULL;
	_arrow = new Arrow(this);
	_vString = Structure_vsh();
	_fString = Structure_fsh();
	this->SlipObject::setName("Segment");
}

void Segment::addMonomerFromAtom(AtomPtr atom)
{
	MonomerPtr mon = atom->getMonomer();
	if (!mon)
	{
		addAtom(atom);
		return;
	}
	
	addAtomsFrom(mon);
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
	clearMesh();
	vec3 c = AtomGroup::centroid();

	for (size_t i = 0; i < atomCount(); i++)
	{
		vec3 a = atom(i)->getAbsolutePosition();
		vec3 diff = vec3_subtract_vec3(c, a);
		vec3_set_length(&diff, 1);
		
		Helen3D::Vertex v;
		memset(&v, '\0', sizeof(Helen3D::Vertex));
		pos_from_vec(&v.pos[0], a);
		pos_from_vec(&v.normal[0], diff);
		
		_vertices.push_back(v);
	}
	
	/* you must now make a mesh and link populate() to other things */
	
	Mesh *m = makeMesh();
	m->setShadersLike(this);
	
	float r = rand() % 360;
	float g = 100;
	float b = 80;
	hsv_to_rgb(r, g, b);

	mesh()->setColour(r, g, b);
	_red = r;
	_green = g;
	_blue = b;
}

void Segment::forceSisterColour()
{
	if (_sister && _sister->hasMesh())
	{
		_sister->mesh()->setColour(_red, _green, _blue);
	}

}

void Segment::refineMesh()
{
	if (_w == NULL)
	{
		_w = new QThread();
	}
	
	mesh()->moveToThread(_w);
	_w->start();
	connect(mesh(), SIGNAL(resultReady()), this, SLOT(handleMesh()));
	connect(this, SIGNAL(refine()), mesh(), SLOT(shrinkWrap()));
	
	emit refine();
}

void Segment::handleMesh()
{
	disconnect(this, SIGNAL(refine()), nullptr, nullptr);
	disconnect(mesh(), SIGNAL(resultReady()), this, SLOT(handleMesh()));
	_w->quit();
	
	mesh()->changeToTriangles();
	
	if (sister() && sister()->ready())
	{
		sister()->_arrow->populate();
		_arrow->populate();
	}
}

bool Segment::ready()
{
	return (!_w || !_w->isRunning());
}

void Segment::render(SlipGL *gl)
{
	glDepthMask(GL_FALSE);
	mesh()->reorderIndices();
	mesh()->render(gl);
	glDepthMask(GL_TRUE);
	_arrow->render(gl);
	glColorMaski(0, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	mesh()->render(gl);
	glColorMaski(0, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}

bool Segment::addAtomIfValid(AtomPtr a)
{
	if (!a)
	{
		return false;
	}

	if (atomCount() == 0)
	{
		addAtom(a);
		return true;
	}

	AtomList cas = findAtoms("CA");
	AtomPtr last = cas[cas.size() - 1];

	if (last->getResidueNum() != a->getResidueNum() - 1)
	{
		return false;
	}

	for (size_t j = 0; j < atomCount(); j++)
	{
		AtomPtr other = atom(j);
		double val = _diff->valueBetween(a, other);
		if (fabs(val) > _threshold || val != val)
		{
			return false;
		}
	}

	addAtom(a);
	
	return true;
}

bool Segment::enoughCommonGround(Segment *other, double mult)
{
	if (!ready() || !other->ready())
	{
		return false;
	}

	double th = _threshold * mult;

	for (size_t j = 0; j < atomCount(); j++)
	{
		AtomPtr aj = atom(j);
		
		for (size_t i = 0; i < other->atomCount(); i++)
		{
			AtomPtr ai = other->atom(i);
			double val = _diff->valueBetween(ai, aj);

			if (fabs(val) > th || val != val)
			{
				return false;
			}
		}
	}
	
	return true;
}

void Segment::addChild(Segment *s)
{
	_children.push_back(s);
	addAtomsFrom(s);
}

void Segment::motionComparedTo(Segment *s, vec3 *start, vec3 *dir)
{
	AtomList cas = findAtoms("CA");
	*start = this->AtomGroup::centroid();
	*dir = empty_vec3();
	double count = 0;

	for (size_t i = 0; i < cas.size(); i++)
	{
		AtomPtr ca = cas[i];
		AtomPtr other = _diff->getPartnerAtom(ca);
		
		if (!other)
		{
			continue;
		}
		
		count++;
		vec3 v1 = ca->getAbsolutePosition();
		vec3 v2 = other->getAbsolutePosition();
		
		vec3 diff = vec3_subtract_vec3(v2, v1);
		vec3_add_to_vec3(dir, diff);
	}

	double mult = 10 / count;
	vec3_mult(dir, mult);
	vec3 half = *dir;
	vec3_add_to_vec3(dir, *start);
	
	vec3_mult(&half, 0.5);
	vec3_subtract_from_vec3(start, half);
	vec3_subtract_from_vec3(dir, half);
}

Segment *Segment::segmentFrom(Segment *s, Segment *t, double mult)
{
	Segment *master = new Segment(s->_diff, s->_threshold * mult);

	master->addChild(s);
	master->addChild(t);
	master->populate();
	
	return master;
}

Segment::~Segment()
{
	for (size_t i = 0; i < _children.size(); i++)
	{
		delete _children[i];
	}

	_children.clear();
	delete _arrow;
}
