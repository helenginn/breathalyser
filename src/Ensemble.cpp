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
#include <algorithm>
#include <libsrc/Polymer.h>
#include <libsrc/Monomer.h>
#include <libsrc/Atom.h>
#include <libinfo/GeomTable.h>
#include <hcsrc/Blast.h>

Ensemble::Ensemble(Ensemble *parent, CrystalPtr c) : QTreeWidgetItem(parent),
SlipObject()
{
	_crystal = c;
	_renderType = GL_LINES;
	_isReference = false;
	findChains();

	initializeOpenGLFunctions();
}

void Ensemble::setName(std::string name)
{
	_name = name;
	updateText();
}

void Ensemble::updateText()
{
	std::string prep;
	prep += (_isReference ? "Ref: " : "");
	prep += _name;
	setText(0, QString::fromStdString(prep));
}

std::string Ensemble::generateSequence(std::string chain)
{
	std::map<int, std::string> resMap;
	AtomList atoms = _crystal->findAtoms("CA", INT_MAX, chain);
	
	int min = INT_MAX;
	int max = -INT_MAX;

	for (size_t i = 0; i < atoms.size(); i++)
	{
		AtomPtr a = atoms[i];
		
		if (!a->getMonomer())
		{
			continue;
		}

		int resi = a->getResidueNum();
		std::string tlc = a->getMonomer()->getIdentifier();
		std::string code = GeomTable::getResCode(tlc); 
		
		resMap[resi] = code;
		
		if (resi < min)
		{
			min = resi;
		}
		
		if (resi > max)
		{
			max = resi;
		}
	}

	std::string seq;
	for (int i = min; i <= max; i++)
	{
		if (resMap.count(i))
		{
			seq += resMap[i];
		}
		else
		{
			seq += ".";
		}
	}
	
	return seq;
}

void Ensemble::findChains()
{
	if (!_crystal)
	{
		return;
	}
	
	_chains.clear();

	for (size_t i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr m = _crystal->molecule(i);
		
		if (!m->isPolymer())
		{
			continue;
		}

		_chains.push_back(m->getChainID());
	}
	
	std::vector<std::string> all = _chains;
	_chains.clear();
	
	for (size_t i = 0; i < all.size(); i++)
	{
		std::string first;
		first += all[i][0];
		
		if (std::find(_chains.begin(), _chains.end(), first) 
		    != _chains.end())
		{
			continue;
		}

		_chains.push_back(first);
		std::string seq = generateSequence(first);
	}
}

void Ensemble::minMaxResidues(std::string ch, int *min, int *max)
{
	*min = INT_MAX;
	*max = -INT_MAX;

	for (size_t i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr m = _crystal->molecule(i);
		
		if (!m->isPolymer())
		{
			continue;
		}

		if (m->getChainID().substr(0, ch.length()) == ch)
		{
			PolymerPtr p = ToPolymerPtr(m);
			if (p->monomerBegin() < *min)
			{
				*min = p->monomerBegin();
			}
			else if (p->monomerEnd() > *max)
			{
				*max = p->monomerEnd();
			}
		}
	}
}

void control_points(vec3 *a, vec3 b, vec3 c, vec3 *d)
{
	vec3 ca = vec3_subtract_vec3(c, *a);
	vec3 bd = vec3_subtract_vec3(b, *d);
	vec3_mult(&ca, 0.4);
	vec3_mult(&bd, 0.4);
	*a = b;
	*d = c;
	vec3_add_to_vec3(a, ca);
	vec3_add_to_vec3(d, bd);
}

vec3 bezier(vec3 p1, vec3 p2, vec3 p3, vec3 p4, double t)
{
	double c1 = (1 - t) * (1 - t) * (1 - t);
	double c2 = 3 * t * (1 - t) * (1 - t);
	double c3 = 3 * t * t * (1 - t);
	double c4 = t * t * t;

	vec3_mult(&p2, c1);
	vec3_mult(&p3, c4);
	vec3_mult(&p1, c2);
	vec3_mult(&p4, c3);
	
	vec3 add = vec3_add_vec3(p1, p2);
	vec3_add_to_vec3(&add, p3);
	vec3_add_to_vec3(&add, p4);
	
	return add;
}

void Ensemble::convertToBezier()
{
	std::vector<Helen3D::Vertex> vs = _vertices;
	std::vector<GLuint> is = _indices;
	_vertices.clear();
	_indices.clear();

	for (int i = 0; i < (int)is.size() - 1; i += 2)
	{
		vec3 p2 = vec_from_pos(vs[is[i]].pos);
		vec3 p3 = vec_from_pos(vs[is[i + 1]].pos);
		
		bool sameBefore = (i <= 1);
		if (!sameBefore)
		{
			sameBefore |= (is[i - 1] != is[i]);
		}

		bool sameAfter = (i >= (int)is.size() - 2);
		
		if (!sameAfter)
		{
			sameAfter = (is[i + 2] != is[i + 1]);
		}

		if (is[i] < 0)
		{
			sameBefore = true;
		}
		
		if (is[i + 2] >= _vertices.size())
		{
			sameAfter = true;
		}
		
		vec3 p1 = p2;
		vec3 p4 = p3;
		
		if (!sameBefore)
		{
			p1 = vec_from_pos(vs[is[i] - 1].pos);
		}
		if (!sameAfter)
		{
			p4 = vec_from_pos(vs[is[i] + 1].pos);
		}
		
		control_points(&p1, p2, p3, &p4);
		
		for (double t = 0; t < 1; t += 0.1)
		{
			if (!sameBefore || t > 0)
			{
				_indices.push_back(_vertices.size());
				_indices.push_back(_vertices.size() + 1);
			}

			vec3 p = bezier(p1, p2, p3, p4, t);
			addVertex(p, NULL);
		}
		
		if (sameAfter)
		{
			_indices.pop_back();
			_indices.pop_back();
		}
	}
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
	for (int i = 0; i < childCount(); i++)
	{
		Ensemble *e = dynamic_cast<Ensemble *>(child(i));
		e->repopulate();
	}
	
	if (!_crystal)
	{
		return;
	}

	_indices.clear();
	_vertices.clear();

	vec3 nanVec = make_vec3(NAN, NAN, NAN);

	for (size_t i = 0; i < _crystal->moleculeCount(); i++)
	{
		MoleculePtr m = _crystal->molecule(i);
		
		if (!m->isPolymer())
		{
			continue;
		}
		
		for (size_t j = 0; j < m->atomCount(); j++)
		{
			AtomPtr a = m->atom(j);
			if (a->getAtomName() == "CA")
			{
				vec3 pos = a->getPDBPosition();
				addCAlpha(pos);
			}
		}

		addCAlpha(nanVec);
	}

	addCAlpha(nanVec);
	
	convertToBezier();

	if (parent() != NULL)
	{
		setColour(0.5, 0.5, 1.);
	}
}

vec3 Ensemble::averagePos()
{
	if (_crystal)
	{
		return centroid();
	}

	vec3 sum = empty_vec3();
	for (int i = 0; i < childCount(); i++)
	{
		Ensemble *e = dynamic_cast<Ensemble *>(child(i));
		vec3 p = e->averagePos();
		vec3_add_to_vec3(&sum, p);
	}
	
	vec3_mult(&sum, 1 / (double)childCount());
	
	if (sum.x != sum.x)
	{
		return empty_vec3();
	}
	
	return sum;
}

void Ensemble::render(SlipGL *gl)
{
	for (int i = 0; i < childCount(); i++)
	{
		Ensemble *e = dynamic_cast<Ensemble *>(child(i));
		e->render(gl);
	}

	SlipObject::render(gl);
}

std::string Ensemble::findMatchingChain(std::string ch, Ensemble *other)
{
	std::string seq = generateSequence(ch);
	int best_mut = INT_MAX;
	std::string best_ch = "";

	for (size_t i = 0; i < other->chainCount(); i++)
	{
		std::string other_ch = other->chain(i);
		int muts, dels;
		std::string otherseq = other->generateSequence(other_ch);
		compare_sequences(seq, otherseq, &muts, &dels);
		
		if (muts == 0 && dels == 0)
		{
			// can't fault that
			return other_ch;
		}
		else if (muts < best_mut)
		{
			best_ch = other_ch;
			best_mut = muts;
		}
	}
	
	if (best_mut > (int)seq.length() / 2)
	{
		std::cout << "Warning: poor identity with best sequence!" << std::endl;
	}
	
	std::cout << "Comparing " << other->name() << " to " 
	<< name() << std::endl;
	std::cout << "Best match for chain " << ch << " is " 
	<< best_ch << std::endl;

	return best_ch;
}
