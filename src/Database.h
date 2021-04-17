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

#ifndef __Database__Database__
#define __Database__Database__

#include <string>
#include <map>
#include <vector>

typedef struct sqlite3 sqlite3;
typedef struct std::map<std::string, std::string> SeqResult;

class Fasta;

class Database
{
public:
	Database(std::string filename);

	int openConnection();
	void closeConnection();
	
	void beginTransaction();
	void endTransaction();

	void importFasta(Fasta *f, bool overwrite = true);
	std::vector<std::string> countryList();

	void query(std::string query);
	
	const std::vector<SeqResult> &results()
	{
		return _results;
	}
private:
	static int callback(void *nu, int argc, char **argv, char **col_names);

	sqlite3 *_db;
	std::string _filename;
	static std::vector<SeqResult> _results;

};

#endif
