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
#include "Segment.h"
#include "Fasta.h"


#include <h3dsrc/shaders/vStructure.h>
#include <h3dsrc/shaders/fStructure.h>
#include <h3dsrc/Icosahedron.h>
#include <h3dsrc/Text.h>
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
	_fastaCount = 0;
	findChains();
	this->SlipObject::setName("Ensemble");
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

vec3 Ensemble::centroidForChain(std::string chain)
{
	std::map<int, std::string> resMap;
	AtomList atoms = _crystal->findAtoms("CA", INT_MAX, chain);
	vec3 sum = empty_vec3();
	
	for (size_t i = 0; i < atoms.size(); i++)
	{
		AtomPtr a = atoms[i];
		vec3 abs = a->getAbsolutePosition();
		
		vec3_add_to_vec3(&sum, abs);
	}
	
	vec3_mult(&sum, 1 / (double)atoms.size());
	
	return sum;
}

std::string Ensemble::generateSequence(std::string chain, int *minRes)
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
			seq += " ";
		}
	}
	
	if (minRes != NULL)
	{
		*minRes = min;
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
	
	for (size_t i = 0; i < segmentCount(); i++)
	{
		segment(i)->render(gl);
	}
	
	for (size_t i = 0; i < _balls.size(); i++)
	{
		_balls[i]->render(gl);
	}
	
	for (size_t i = 0; i < _texts.size(); i++)
	{
		_texts[i]->render(gl);
	}

	SlipObject::render(gl);
}

std::string Ensemble::findMatchingChain(std::string ch, Ensemble *other)
{
	std::string seq = generateSequence(ch);
	int best_mut = INT_MAX;
	double best_length = FLT_MAX;
	std::string best_ch = "";
	vec3 this_c = centroidForChain(ch);

	for (size_t i = 0; i < other->chainCount(); i++)
	{
		std::string other_ch = other->chain(i);
		int muts, dels;
		std::string otherseq = other->generateSequence(other_ch);
		compare_sequences(seq, otherseq, &muts, &dels);
		std::cout << "Potential sequence: " << muts << " mutations." << std::endl;
		vec3 diff = other->centroidForChain(other_ch);
		vec3_subtract_from_vec3(&diff, this_c);
		double l = vec3_length(diff);
		std::cout << "Potential sequence: " << l << " Å separation." << std::endl;
		
		if (muts < best_mut)
		{
			best_ch = other_ch;
			best_mut = muts;
			best_length = l;
		}
		else if (muts <= (int)seq.length() / 2)
		{
			vec3 diff = other->centroidForChain(other_ch);
			vec3_subtract_from_vec3(&diff, this_c);
			double l = vec3_length(diff);

			if (l < best_length)
			{
				best_ch = other_ch;
				best_mut = muts;
				best_length = l;
			}

		}
	}
	
	if (best_mut > (int)seq.length() / 2)
	{
		std::cout << "Warning: poor identity with best sequence!" << std::endl;
	}
	std::cout << best_mut << " mutations." << std::endl;
	std::cout << best_length << " Å separation." << std::endl;
	
	std::cout << "Comparing " << other->name() << " to " 
	<< name() << std::endl;
	std::cout << "Best match for chain " << ch << " is " 
	<< best_ch << std::endl;

	return best_ch;
}

void Ensemble::deleteSegments()
{
	for (size_t i = 0; i < segmentCount(); i++)
	{
		delete segment(i);
	}
	
	_segments.clear();
}

void Ensemble::processNucleotides(Fasta *f)
{
	int minRes = 0; int maxRes = 0;
	minMaxResidues(chain(0), &minRes, &maxRes);
	f->setOffset(minRes);

	if (f->hasResult())
	{
		return;
	}

	std::string seq = generateSequence(chain(0));
	f->roughCompare(seq, minRes);
}

bool Ensemble::processFasta(Fasta *f, std::string requirements)
{
	bool should = true;
	
	if (f->isProblematic())
	{
		return false;
	}
	
	if (requirements.length() > 0)
	{
		std::vector<std::string> reqs = split(requirements, ',');
		
		for (size_t j = 0; j < reqs.size(); j++)
		{
			bool found = false;
			bool invert = false;
			std::string req = reqs[j];
			
			if (req[0] == '!')
			{
				req.erase(req.begin());
				invert = true;
			}
			
			if (req[0] < '0' || req[0] > '9')
			{
				req.erase(req.begin());
			}
			
			for (size_t i = 0; i < f->mutationCount(); i++)
			{
				std::string variant = f->mutation(i);
				std::string ending = variant.substr(1);
				
				int comp_length = std::min(ending.length(), req.length());

				if (req.substr(0, comp_length) == 
				    ending.substr(0, comp_length))
				{
					found = true;
					break;
				}
			}

			if ((!found && !invert) || (found && invert))
			{
				should = false;
			}
		}
	}
	
	if (!should)
	{
		return false;
	}
	
	for (size_t i = 0; i < f->mutationCount(); i++)
	{
		processMutation(f->mutation(i));
	}
	
	_fastaCount++;

	return true;
}

void Ensemble::processMutation(std::string mutation)
{
	int mut = atoi(&mutation.c_str()[1]);
	AtomList atoms = _crystal->findAtoms("CA", mut);
	
	if (atoms.size() == 0)
	{
		return;
	}
	
	_muts[mut].push_back(mutation);
}

size_t Ensemble::makeBalls()
{
	AtomList atoms = _crystal->findAtoms("CA");
	
	if (_fastaCount <= 1)
	{
		return 0;
	}

	for (size_t i = 0; i < atoms.size(); i++)
	{
		int resNum = atoms[i]->getResidueNum();
		
		if (_muts[resNum].size() == 0)
		{
			continue;
		}
		
		std::string aas;
		unsigned char refaa = _muts[resNum][0][0];
		for (size_t j = 0; j < _muts[resNum].size(); j++)
		{
			unsigned char back = _muts[resNum][j].back();
			
			if (aas.find(back) == std::string::npos)
			{
				aas.push_back(back);
			}
		}

		double counts = _muts[resNum].size();

		double pct = 100 * counts / (double)_fastaCount;
		double pct100 = 10 * pct;
		
		double inflate = log(pct100) / log(10);
		
		/* no negatives! */
		if (inflate < 0)
		{
			continue;
		}

		vec3 abs = atoms[i]->getAbsolutePosition();

		if (pct > 0.3)
		{
			std::string str;
			str += refaa;
			str += i_to_str(resNum) + aas + ", ";
			str += f_to_str(pct, 1) + "%";
			Text *text = new Text();
			text->setProperties(abs, str, 64, Qt::black,
			                    0, 20, -20);
			text->prepare();
			_texts.push_back(text);
			_textMap[atoms[i]] = text;
		}

		Icosahedron *ico = new Icosahedron();
		ico->setPosition(abs);
		ico->setColour(1, 0.3, 0.3);
		
		if (aas[0] == '-')
		{
			ico->setColour(0.3, 0.3, 0.3);
		}
		else if (aas[0] == '+')
		{
			ico->setColour(0.3, 0.3, 1.0);
		}

		ico->triangulate();
		std::string v = Structure_vsh();
		std::string f = Structure_fsh();
		ico->changeProgram(v, f);
		ico->resize(inflate);

		_balls.push_back(ico);
		_ballMap[atoms[i]] = ico;
	}
	
	return _fastaCount;
}

void Ensemble::clearBalls()
{
	for (size_t i = 0; i < _balls.size(); i++)
	{
		delete _balls[i];
	}

	for (size_t i = 0; i < _texts.size(); i++)
	{
		delete _texts[i];
	}

	_texts.clear();
	_balls.clear();
	_ballMap.clear();
	_textMap.clear();
	_muts.clear();
	_fastaCount = 0;
}
