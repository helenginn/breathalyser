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

#ifndef __breathalyser__ensemble__ 
#define __breathalyser__ensemble__ 

#include <QTreeWidgetItem>
#include <h3dsrc/SlipObject.h>

class Text;
class Fasta;
class Segment;
class Icosahedron;

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;

/******************** Crystal definition ***************************/
namespace Vagabond
{
	class Crystal;
}

#define ToCrystalPtr(a) (boost::static_pointer_cast<Vagabond::Crystal>((a)))
typedef boost::shared_ptr<Vagabond::Crystal> CrystalPtr;
typedef boost::weak_ptr<Vagabond::Crystal> CrystalWkr;

/****************** end Crystal definition *************************/

class Ensemble : public QTreeWidgetItem, public SlipObject
{
public:
	Ensemble(Ensemble *parent, CrystalPtr c);
	
	void setName(std::string name);
	
	std::string name()
	{
		return _name;
	}
	
	size_t chainCount()
	{
		return _chains.size();
	}
	
	std::string chain(int i)
	{
		return _chains[i];
	}

	CrystalPtr crystal()
	{
		return _crystal;
	}
	
	void setReference(bool ref)
	{
		_isReference = ref;
		updateText();
	}
	
	void addSegment(Segment *seg)
	{
		_segments.push_back(seg);
	}
	
	size_t segmentCount()
	{
		return _segments.size();
	}
	
	Segment *segment(int i)
	{
		return _segments[i];
	}
	
	void setSegments(std::vector<Segment *> segments)
	{
		_segments = segments;
	}
	
	void removeSegment(int i)
	{
		_segments.erase(_segments.begin() + i);
	}
	
	std::vector<Segment *> segments()
	{
		return _segments;
	}
	
	void deleteSegments();

	void processNucleotides(Fasta *f);
	void addCAlpha(vec3 point);
	void repopulate();
	void updateText();
	vec3 averagePos();

	virtual void render(SlipGL *gl);

	size_t makeBalls();
	bool processFasta(Fasta *f, std::string requirements = "");
	void clearBalls();

	void minMaxResidues(std::string ch, int *min, int *max);
	std::string generateSequence(std::string chain, int *minRes = NULL);
	vec3 centroidForChain(std::string chain);
	std::string findMatchingChain(std::string ch, Ensemble *other);
	
private:
	void processMutation(std::string mutation);
	std::vector<Segment *> _segments;
	std::vector<Icosahedron *> _balls;
	std::vector<Text *> _texts;
	std::map<AtomPtr, Icosahedron *> _ballMap;
	std::map<AtomPtr, Text *> _textMap;
	std::map<int, std::vector<std::string> > _muts;
	int _fastaCount;
	void findChains();
	void convertToBezier();
	CrystalPtr _crystal;

	bool _isReference;
	std::string _name;
	std::vector<std::string> _chains;
};

#endif
