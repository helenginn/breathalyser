// splitseq
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

#ifndef __splitseq__Fetch__
#define __splitseq__Fetch__

#include <QMainWindow>
#include "Database.h"

class Fetch : public QMainWindow
{
Q_OBJECT
public:
	Fetch(QWidget *parent);

	void setDatabase(Database *db)
	{
		_db = db;
		populate();
	}

	void populate();

private slots:
	void getSequences();
	void updateCount();
private:
	void processResults(std::vector<SeqResult> results);
	std::string makeQuery(bool how_many);

	Database *_db;
};

#endif
