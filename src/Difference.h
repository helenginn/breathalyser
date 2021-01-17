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

#ifndef __breathalyser__difference__
#define __breathalyser__difference__

#include <QImage>
#include <QTreeWidgetItem>
#include <hcsrc/vec3.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class StructureView;
class DiffDisplay;
class Ensemble;
class Segment;
class Main;

/****************** Atom definition ********************/

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;
typedef boost::weak_ptr<Atom> AtomWkr;

/****************** End atom definition ********************/

typedef std::pair<AtomPtr, AtomPtr> AtomCouple;
typedef std::pair<int, int> IntPair;

class Difference : public QObject, public QTreeWidgetItem, public QImage
{
Q_OBJECT
public:
	Difference(int w, int h);

	void toCoupleView(StructureView *view);
	void setEnsembles(Ensemble *a, Ensemble *b);
	void setMain(Main *m)
	{
		_main = m;
	}
	
	void setDisplay(DiffDisplay *display)
	{
		_display = display;
	}
	
	size_t atomCount()
	{
		return _atoms.size();
	}
	
	size_t segmentCount()
	{
		return _segments.size();
	}
	
	Segment *segment(int i)
	{
		return _segments[i];
	}

	void populate(bool force = false);
public slots:
	void thresholdChanged(int val);
	void calculate();

private:
	void findCommonAtoms();
	void findSegments(double val);
	void loadSegmentsToView();

	Ensemble *_ea;
	Ensemble *_eb;
	
	Main *_main;

	DiffDisplay *_display;
	StructureView *_coupleView;
	std::vector<AtomPtr> _atoms;
	std::vector<AtomCouple> _atomCouples;
	std::vector<Segment *> _segments;
	std::map<AtomCouple, double> _vals;
	double _max;
	bool _drawn;
};

#endif
