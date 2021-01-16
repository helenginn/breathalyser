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

class Ensemble;

/****************** Atom definition ********************/

class Atom;
typedef boost::shared_ptr<Atom> AtomPtr;
typedef boost::weak_ptr<Atom> AtomWkr;

/****************** End atom definition ********************/

typedef std::pair<AtomPtr, AtomPtr> AtomCouple;

class Difference : public QTreeWidgetItem, public QImage
{
public:
	Difference(int w, int h);

	void setEnsembles(Ensemble *a, Ensemble *b);

	void populate();
private:
	void findCommonAtoms();

	Ensemble *_ea;
	Ensemble *_eb;
	
	std::vector<AtomCouple> _atomCouples;
};

#endif
