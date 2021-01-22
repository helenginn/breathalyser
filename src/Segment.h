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

#ifndef __breathalyser__segment__
#define __breathalyser__segment__

#include <libsrc/AtomGroup.h>
#include <h3dsrc/SlipObject.h>
#include <QObject>

class Arrow;
class QThread;
class Difference;

class Segment : public QObject, public SlipObject, public AtomGroup
{
Q_OBJECT
public:
	Segment(Difference *d, double threshold);
	static Segment *segmentFrom(Segment *s, Segment *t, double mult);
	virtual ~Segment();
	
	bool ready();
	
	void startEnd(int *min, int *max);
	void populate();

	virtual void render(SlipGL *gl);
	
	void addMonomerFromAtom(AtomPtr);
	void kickOutEarly();
	bool addAtomIfValid(AtomPtr a);
	void forceSisterColour();
	bool enoughCommonGround(Segment *other, double mult);

	size_t subCount()
	{
		if (_children.size() == 0)
		{
			return 1;
		}

		return _children.size();
	}
	
	Segment *sub(int i)
	{
		if (_children.size() == 0)
		{
			return this;
		}
		else return _children[i];
	}
	
	void setSister(Segment *s)
	{
		_sister = s;
	}
	
	Segment *sister()
	{
		return _sister;
	}

	void motionComparedTo(Segment *s, vec3 *start, vec3 *dir);
signals:
	void refine();
public slots:
	void refineMesh();
	void handleMesh();

protected:

private:
	void addChild(Segment *s);

	Arrow *_arrow;
	QThread *_w;
	Segment *_sister;
	Difference *_diff;
	double _threshold;
	std::vector<Segment *> _children;
};

#endif
