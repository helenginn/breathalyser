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

#include "Ensemble.h"
#include <iostream>
#include <c4xsrc/QuickAtoms.h>
#include "CTrace.h"

Ensemble::Ensemble(Ensemble *parent, QuickAtoms *qa) : QTreeWidgetItem(parent)
{
	_qa = qa;
	_renderType = GL_LINES;
}

void Ensemble::setName(std::string name)
{
	_name = name;
	setText(0, QString::fromStdString(_name));
}

void control_points(vec3 *a, vec3 b, vec3 c, vec3 *d)
{
	vec3 ba = vec3_subtract_vec3(b, *a);
	vec3 dc = vec3_subtract_vec3(*d, c);
	vec3 bc = vec3_subtract_vec3(b, c);
	vec3 cb = vec3_subtract_vec3(c, b);
	
	vec3 abc = vec3_add_vec3(ba, cb);
	vec3 dcb = vec3_add_vec3(dc, bc);
	vec3_mult(&abc, 0.67);
	vec3_mult(&dcb, 0.67);
	vec3_add_to_vec3(a, abc);
	vec3_add_to_vec3(d, dcb);
}

void Ensemble::addCAlpha(vec3 point)
{
	bool lastOK = (_vertices.size() > 0);
	
	if (lastOK)
	{
		Helen3D::Vertex last = _vertices.back();
		lastOK = (last.pos[0] == last.pos[0]);
	}
	
	if (point.x == point.x && lastOK)
	{
		_indices.push_back(_vertices.size() - 1);
		_indices.push_back(_vertices.size());
	}

	Helen3D::Vertex v;
	memset(v.pos, 0, sizeof(GLfloat) * 3);
	
	v.color[3] = 1;
	v.pos[0] = point.x;
	v.pos[1] = point.y;
	v.pos[2] = point.z;

	_vertices.push_back(v);

}

void Ensemble::repopulate()
{
	_indices.clear();
	_vertices.clear();

	vec3 nanVec = make_vec3(NAN, NAN, NAN);

	size_t count = _qa->chainCount();
	std::cout << "Count: " << count << std::endl;

	for (size_t i = 0; i < count; i++)
	{
		std::string chain = _qa->chain(i);
		Vec3Vec poz = _qa->posForChain(chain);
		
		for (size_t j = 0; j < poz.size(); j++)
		{
			vec3 pos = poz[j];
			addCAlpha(pos);
		}

		addCAlpha(nanVec);
	}

	addCAlpha(nanVec);
}
