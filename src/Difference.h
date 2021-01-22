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

class CoupleDisplay;
class DiffDisplay;
class Ensemble;
class Segment;
class Main;

/****************** Atom definition ********************/

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;
typedef boost::weak_ptr<Atom> AtomWkr;

/****************** End atom definition ********************/

typedef std::map<AtomPtr, AtomPtr> AtomMap;
typedef std::pair<AtomPtr, AtomPtr> AtomCouple;
typedef std::pair<int, int> IntPair;

class Difference : public QObject, public QTreeWidgetItem, public QImage
{
Q_OBJECT
public:
	Difference(int w, int h);

	void toCoupleDisplay(CoupleDisplay *view);
	void setEnsembles(Ensemble *a, Ensemble *b);
	void applySegmentsToEnsembles();

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
	
	Ensemble *aEnsemble()
	{
		return _ea;
	}
	
	Ensemble *bEnsemble()
	{
		return _eb;
	}

	void populate(bool force = false);
	double valueBetween(AtomPtr a, AtomPtr b);
	
	AtomPtr getPartnerAtom(AtomPtr a);
public slots:
	void thresholdChanged(int val);
	void calculate();
	void tryMerges();

private:
	void drawOnSegment(Segment *s, Segment *t, 
	                   QPainter &painter, double box_size);

	void findCommonAtoms();
	void findSegments(double val);
	void addSegments(Segment *seg_a, double threshold);
	bool closerSegmentThan(Segment *s, Segment *t);

	Ensemble *_ea;
	Ensemble *_eb;

	std::vector<Segment *> _aSegments;
	std::vector<Segment *> _bSegments;
	
	Main *_main;

	DiffDisplay *_display;
	CoupleDisplay *_coupleDisplay;
	std::vector<AtomPtr> _atoms;
	std::vector<AtomCouple> _atomCouples;
	AtomMap _map;

	std::map<AtomCouple, double> _vals;
	double _max;
	bool _drawn;
};

#endif
