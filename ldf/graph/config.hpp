/*
 * Configuration.hpp
 *
 *  Created on: Nov 05, 2010
 *      Author: Thang N. Dinh
 *      Version: 1.0 Nov 05, 2010
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
#include <limits>

#ifndef     CONFIGURATION_H_
#define     CONFIGURATION_H_

/**
 * @brief Configuration loader
 * The syntax for the Configuration file is very simple
 * 1) Comment start with '#'
 * 2) Each line has syntax:   <key>=<value>
 * 3) Keys are insensitive and all space characters are ignored. For example, "number of vertices",
 * "numberofvertices","Number of vertices" are all treated as  a same key.
 * 4) One key can associate with many values.
 *
 * REMARK: When a key is associated with multiple values, the starting index of the values 
 * is ONE (not ZERO). For example, to access the second value that associates with key
 * "input", you invoke the function value("input",2).
 */
class Configuration {
	std::vector<std::string> keys, values;
	std::istringstream is;
public:
	/*
	 * @brief Load options from file
	 */
	bool load(std::string fileName) {
		std::ifstream fi(fileName.c_str());
		if (!fi.good()) {
			std::cerr << "Error: Loading configuration from file " << fileName
					<< " failed!" << std::endl;
			return false;
		}
		std::string s;
		int line = 0;
		while (!fi.eof()) {
			line++;
#undef max
			if (fi.peek() == '#' || fi.eof()) {
				fi.ignore(std::numeric_limits<int>::max(), '\n');
				continue;
			}
			getline(fi, s);
			size_t eq = s.find_first_of('=');
			if (eq == std::string::npos) {
				if (s.find_first_not_of(' ') != std::string::npos)
					std::cerr << "Warning: Invalid setting at line " << line
							<< ": " << s << std::endl;
				continue;
			}
			std::string key = s.substr(0, eq), value = s.substr(eq + 1);
			keys.push_back(normalized(key));
			values.push_back(value);
		} // while
		fi.close();
		return true;
	}


	/**
	 *
	 * @param s	input string
	 * @return A string without space characters
	 */
	static std::string normalized(const std::string &s) {
		// To lower case
		std::string result = "";
		for (size_t i = 0; i < s.length(); ++i)
			if (s[i] != ' ' && s[i] != '\t' && s[i] != '\n')
				result += tolower(s[i]);
		return result;
	}

	/**
	 * @param key	A key to look up
	 * @return The value associated with the key or a stream of emtpy if the key doesn't exist
	 */
	std::istringstream &operator[](std::string key) {
		key = normalized(key);
		for (size_t i = 0; i < keys.size(); ++i)
			if (key == keys[i]) {
				is.clear(); // Reset the state to avoid EOF
				is.str(values[i]);
				//              cout <<"Key "<<key<<"'s value: "<< is.str() << endl;
				return is;
			}
		std::cerr << "Warning: Key " << key << " not found!" << std::endl;
		is.str(" ");
		return is;
	}
	/**
	 * Retrieve the value associated with a key when a key has multiple values.
	 * @param key		The key
	 * @param order		The order in which the key appear. The order starts from 1.
	 */
	std::string value(std::string key, int order = 1) {
		key = normalized(key);
		int count = 0;
		for (size_t i = 0; i < keys.size(); ++i)
			if (key == keys[i]) {
				count++;
				if (count == order)
					return values[i];
			}
		std::cerr << "Warning: Key " << key << ", order: " << order
				<< " not found!" << std::endl;
		return "";
	}

	/**
	 * Retrieve the value associated with a key when a key has multiple values.
	 * @param key		The key
	 * @param order		The order in which the key appear. The order starts from 1.
	 */
	template <class T>
	T cast(std::string key, int order = 1) {
		T x;
		std::istringstream is(value(key, order));
		is >> x;
		return x;
	}

	/**
	 * Compare the value associated with a key with a given value.
	 * Both strings are normalized before comparison.
	 * @param key		The key
	 * @param value		The string to compare to.
	 * @param order 	The order in which key appears (multiple values)
	 * @return	true	if two values are equal after normalization
	 * 			false	if key does not match or exist
	 */
	bool match(const std::string &key, const std::string &vl, int order = 1) {
		std::string s1 = normalized(value(key, order)), s2 = normalized(vl);
		return s1 == s2;
	}

	/**
	 *
	 * @param key	Key to look up
	 * @return The number of values associate with the given key
	 */
	int count(std::string key) const  {
		key = normalized(key);
		int count = 0;
		for (size_t i = 0; i < keys.size(); ++i)
			if (key == keys[i])
				count++;
		return count;
	}
	/**
	 * Retrieve the value associated with a key in form of a stream
	 * @param key		The key
	 * @param order		The order in which the key appear. The order starts from 1.
	 */
	std::istringstream &valueStream(std::string key, int order = 1) {
		is.clear(); // Reset the stream's state to clear EOF
		is.str(value(key, order));
		return is;
	}
	/**
	 *
	 * @param s	 The input string
	 * @return A string obtained by removal of all space characters in s
	 */
	static std::string trim(const std::string &s) {
		int start = s.find_first_not_of(" \t\n"), len = s.find_last_not_of(
				" \t\n") - start + 1;
		return s.substr(start, len);
	}
	/**
	 * @brief Check if a key exists
	 */
	bool exists(const std::string &key) const {
		return this->count(key) > 0;
	}

	/**
	 * @brief Add a pair of (key, value) to the configuration
	 * @param key
	 * @param value
	 */
	void add(const std::string &key, const std::string &value)  {
		keys.push_back(normalized(key));
		values.push_back(value);
	}

	/**
	 * @brief Load the Configuration from a list of (program's) parameters
	 * @param argc	The number of arguments
	 * @param args  Array of arguments
	 */
	Configuration(int argc, char* args[]) {
		int p = 1;
		while (p < argc) {
			if (args[p][0] != '-') {
				std::cerr
						<< "Error: Invalid syntax. Please, see help! Syntax: "
						<< "<program> -help" << std::endl;
				exit(1);
			}
			string key = args[p];
			key = key.substr(1);
			string value = "";
			for (++p; (p < argc) &&(args[p][0] != '-'); ++p) {
					if (value != "")
						value += string(" ");
					value += string(args[p]);
		    }
			keys.push_back(normalized(key));
			values.push_back(value);
		}
		int cfg_file = this->count("config");
		// Load all provided Configuration file(s)
		if (cfg_file > 0) {
			for (int i = 1; i <= cfg_file; ++i)
				this->load(value("config", i));
		}
	}

	/**
	 * Create an alternative name for a key.
	 * @param key The key's name
	 * @param alias The alias name of the key
	 * @return "true" if the key exists, "false" otherwise
	 */
	bool alias(const std::string &key, const std::string alias) {
		int cnt = count(key);
		for(int i=1; i <= cnt;++i)
			add( alias, value(key, i));
		return (cnt >0);
	}

	/**
	 * Load the Configuration from a given file
	 * @param fileName	The Configuration file
	 */
	Configuration(std::string fileName) {
		load(fileName);
	}

};

#endif

//int main() {
//	Configuration cf("test.txt");
//	std::string s1, s2;
//	cf["name"] >> s1;
//	std::cout << "'" << s1 << "'" << std::endl;
//	std::cout << "'" << cf["life"].str() << "'" << std::endl;
//	for (int i = 0; i < cf.count("na me"); ++i)
//		std::cout << i << "\t" << cf.value("name", i + 1) << std::endl;
//	std::cout << "'" << cf.trim(" asaf asfd ewr ewr                       ") << "'"
//			<< std::endl;
//	system("pause");
//}
