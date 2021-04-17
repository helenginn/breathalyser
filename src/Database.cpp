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

#include "Database.h"
#include "Fasta.h"
#include <sqlite3.h>
#include <hcsrc/FileReader.h>

std::vector<SeqResult> Database::_results;

Database::Database(std::string filename)
{
	_db = NULL;
	_filename = filename;
}

int Database::openConnection()
{
//	std::cout << _filename << std::endl;
	int rc = sqlite3_open(_filename.c_str(), &_db);

	if (rc)
	{
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(_db));
		sqlite3_close(_db);
		return(1);
	}

	return 0;
}

int Database::callback(void *nu, int argc, char **argv, char **col_names)
{
	SeqResult result;

	for (int i = 0; i < argc; i++)
	{
		result[col_names[i]] = (argv[i] != NULL ? argv[i] : "NULL");
	}

	_results.push_back(result);
	return 0;
}

void Database::query(std::string query)
{
	char *zErrMsg = 0;
	_results.clear();

	int rc = sqlite3_exec(_db, query.c_str(), callback, 0, &zErrMsg);

	if (rc != SQLITE_OK)
	{
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	}
}

void Database::closeConnection()
{
	if (_db != NULL)
	{
		sqlite3_close(_db);
		_db = NULL;
	}
}

void Database::importFasta(Fasta *f, bool overwrite)
{
//	std::string q = f->selectQuery();
	std::string q = f->insertQuery();
	q += f->updateQuery();
	query(q);
	return;
	
	if (_results.size() == 0)
	{
		std::string q = f->insertQuery();
		query(q);
		std::cout << "Inserted " << f->name() << std::endl;
	}
	else if (overwrite)
	{
		std::string q = f->updateQuery();
		query(q);
		std::cout << "Updated " << f->name() << std::endl;
	}
	else
	{
		std::cout << "Skipping " << f->name() << std::endl;
	}
}


std::vector<std::string> Database::countryList()
{
	openConnection();

	std::string q = "SELECT DISTINCT country FROM sequences ORDER BY country;";
	query(q);

	std::vector<std::string> list;
	for (size_t i = 0; i < _results.size(); i++)
	{
		list.push_back(_results[i]["country"]);
	}
	
	closeConnection();
	return list;
}

void Database::beginTransaction()
{
	std::string q = "BEGIN TRANSACTION;";
	query(q);
}

void Database::endTransaction()
{
	std::string q = "END TRANSACTION;";
	query(q);
}
