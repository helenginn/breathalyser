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

#include "DiffDisplay.h"
#include "Difference.h"

DiffDisplay::DiffDisplay(QWidget *parent, Difference *diff) : QLabel(parent)
{
	_diff = diff;
}

void DiffDisplay::changeDifference(Difference *diff)
{
	if (diff != NULL)
	{
		_diff = diff;
	}
	
	if (_diff == NULL)
	{
		return;
	}

	QImage i = diff->scaled(width(), height(), Qt::IgnoreAspectRatio);
	setPixmap(QPixmap::fromImage(i));

}
