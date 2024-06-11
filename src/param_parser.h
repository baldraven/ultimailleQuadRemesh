#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <typeinfo>

// https://github.com/ultimaille/param-parser

struct Parameters {

	// List of default param types (like an "enum")
	struct Type {
	#define SCALAR_PARAM_TYPE(x, s) static constexpr const char* x {s};
	#define VECTOR_PARAM_TYPE(x, s) static std::string x(int dim) { return std::string(s) + "." + std::to_string(dim); }

		// Primitive types
		SCALAR_PARAM_TYPE(Int, "int")
		SCALAR_PARAM_TYPE(Float, "float")
		SCALAR_PARAM_TYPE(Double, "double")
		SCALAR_PARAM_TYPE(Bool, "bool")
		SCALAR_PARAM_TYPE(String, "string")
		SCALAR_PARAM_TYPE(File, "file")
		SCALAR_PARAM_TYPE(Enum, "enum")
		SCALAR_PARAM_TYPE(Input, "input")

		// Mesh types
		SCALAR_PARAM_TYPE(Mesh, "mesh")
		SCALAR_PARAM_TYPE(MeshTri, "mesh.tri")
		SCALAR_PARAM_TYPE(MeshQuad, "mesh.quad")
		SCALAR_PARAM_TYPE(MeshTet, "mesh.tet")
		SCALAR_PARAM_TYPE(MeshHex, "mesh.hex")

		// Attribute types
		VECTOR_PARAM_TYPE(VerticesInt, "vertices.int")
		VECTOR_PARAM_TYPE(VerticesFloat, "vertices.float")
		VECTOR_PARAM_TYPE(VerticesDouble, "vertices.double")
		VECTOR_PARAM_TYPE(VerticesBool, "vertices.bool")

		VECTOR_PARAM_TYPE(FacetsInt, "facets.int")
		VECTOR_PARAM_TYPE(FacetsFloat, "facets.float")
		VECTOR_PARAM_TYPE(FacetsDouble, "facets.double")
		VECTOR_PARAM_TYPE(FacetsBool, "facets.bool")

		VECTOR_PARAM_TYPE(EdgesInt, "edges.int")
		VECTOR_PARAM_TYPE(EdgesFloat, "edges.float")
		VECTOR_PARAM_TYPE(EdgesDouble, "edges.double")
		VECTOR_PARAM_TYPE(EdgesBool, "edges.bool")

		VECTOR_PARAM_TYPE(CellsInt, "cells.int")
		VECTOR_PARAM_TYPE(CellsFloat, "cells.float")
		VECTOR_PARAM_TYPE(CellsDouble, "cells.double")
		VECTOR_PARAM_TYPE(CellsBool, "cells.bool")
	
	#undef SCALAR_PARAM_TYPE
	#undef VECTOR_PARAM_TYPE

	};

	// List of default kind of params (like an "enum")
	struct Kind {
		static constexpr const char* basic {"basic"};
		static constexpr const char* advanced {"advanced"};
		static constexpr const char* system {"system"};
	};

#define ADD_FIELD(x) Param& x(std::string str) {  _##x=str; return *this; }\
	std::string _##x;
	
	struct Param {

		Param() {
			_type = "undefined";
			_value = "undefined";
			_description = "undefined";
			_possible_values = "undefined";
			_type_of_param = "basic";
			_visible = "true";
		}

		ADD_FIELD(type)				// should be one of "int","float","string","file","directory","bool","enum"
		ADD_FIELD(value)
		ADD_FIELD(description)
		ADD_FIELD(possible_values)		// string with values are separated by ','. For example "option1,option2,option3"
		ADD_FIELD(type_of_param)		// should be one of "basic", "advanced" or "system"
		ADD_FIELD(visible)		// should be true or false
		
		Param & default_value(std::string str) {  _value=str; return *this; }

		
		operator int() {
			assert_type_equals(Type::Int);
			try { return std::stoi(_value); }
			catch (std::invalid_argument&) { throw_type_cast_error(_value, Type::Int); return -1; }
		}
		
		operator bool() {
			assert_type_equals(Type::Bool);
			return is("true");
		}

		operator float() {
			assert_type_equals(Type::Float);
			try { return std::stof(_value); }
			catch (std::invalid_argument&) { throw_type_cast_error(_value, Type::Float); return -1; }
		}
		operator double() {
			assert_type_equals(Type::Double);
			try { return std::stod(_value); }
			catch (std::invalid_argument&) { throw_type_cast_error(_value, Type::Double); return -1; }
		}
		operator std::string() { return _value; }

		
		inline bool is(std::string str) { return _value.compare(str) == 0; }

		// Eventually format value (e.g: quotes for string)
		inline std::string formatted_value() {
			if (_type.compare(Type::String) == 0)
				return "\"" + _value + "\"";
			else 
				return _value;
		}

		// Set value
		inline void set(std::string value) {
		    _value = value;
		}

		private:
			inline void assert_type_equals(const std::string type) {
				if (_type.compare(type) != 0)
					throw_type_cast_error(_type, type);
			}

			inline void throw_type_cast_error(const std::string source_type, const std::string target_type) {
				throw std::runtime_error("Unable to cast " + source_type + " into an " + target_type + ".");
			}

	};
#undef ADD_FIELD

	std::string help;

    // Default constructor
    Parameters() {
		// Declare automatically special parameters
		this->add("string", "result_path", "").type_of_param(Parameters::Kind::system);
		this->add("string", "run_from", "").type_of_param(Parameters::Kind::system);
	}

    // Init with serialized parameters string
    Parameters(std::string serialized_params) {

        std::istringstream s(serialized_params);
        std::string line;

        const std::string delimiter = ";";

        while (std::getline(s, line)) {

            std::size_t last_pos = 0;
            // Assume there is only 6 features per parameter
            // Doesn't check whether file is well formed
            std::string chunks[6];

            // Skip comments (begin by char #)
            if (line[0] == '#')
                continue;

            std::size_t pos = 0;
            std::map<std::string, std::string> kv;
            while ((pos = line.find(delimiter, last_pos)) != std::string::npos) {
                std::string chunk = line.substr(last_pos, pos - last_pos);
                
                std::size_t equal_pos = chunk.find("=");
                if (equal_pos == std::string::npos)
                    throw std::runtime_error("Malformed parameter string.");

                std::string key = chunk.substr(0, equal_pos);
                std::string val = chunk.substr(equal_pos + 1, chunk.length() - equal_pos - 1);

                kv[key] = val;
                last_pos = pos + 1;
            }

            // Add param
            auto & param = add(kv["type"], kv["name"], kv["value"]);
            param._possible_values = kv["possible_values"];
            param._description = kv["description"];
            param._type_of_param = kv["type_of_param"];
        }

    }

	Param& add(std::string type, std::string name, std::string default_value) {
		// if (data.find(name) != data.end())
		// 	throw std::runtime_error("Duplicate parameter '" + name + "' found.");

		data[name] = Param();
		data[name]._type = type;
		data[name]._value = default_value;

		return data[name];
	}

	Param& operator[](std::string name) { return data[name]; }

	void init_from_string(const std::string& s) {
		int eq_pos = 0;
		for (std::size_t i = 0; i < s.size();i++) if (s[i] == '=') eq_pos = i;
		std::string var_name = s.substr(0, eq_pos);
		if (data.find(var_name) == data.end())
			std::cerr << "Argument " << var_name << " is not defined for this binary\n";
		else {
			if (s.substr(eq_pos + 1, s.size()).compare("nil")!=0) 
				data[s.substr(0, eq_pos)]._value = s.substr(eq_pos + 1, s.size());
		}
	}

    std::string str_values() {
        std::string s;
		for (auto it : data)
			s += it.first + "=" + it.second.formatted_value() + " ";

        return s;
    }

	void show_values() {
	    std::cerr << str_values();
	}

	void init_from_args(int argc, char** argv) {
		// export parameters "API"
		if (argc > 1 && std::string("--show-params").compare(argv[1]) == 0) {
			std::cout << "#This file contains reflexion information for calling a binary file\n";
			for (auto it : data) {
			std::cout
			    << "name=" << it.first << ";"
			    << "type=" << it.second._type << ";"
			    << "value=" << it.second._value << ";"
			    << "possible_values=" << it.second._possible_values << ";"
			    << "description=" << it.second._description << ";"
			    << "type_of_param=" << it.second._type_of_param << ";"
			    << "visible=" << it.second._visible
			    << std::endl;
			}
			exit(EXIT_SUCCESS);
		}
		else if (argc > 1 && std::string("-h").compare(argv[1]) == 0) {
			std::cout << help << std::endl;
			exit(EXIT_SUCCESS);
		}

		// read arguments from command line
		for (int p = 0; p < argc - 1;p++) init_from_string(argv[p + 1]);
	}

	// Special parameters
	std::string result_path() {
		return data["result_path"]._value;
	}

	std::string run_from() {
		return data["run_from"]._value;
	}

	bool has_result_path() {
		return !result_path().empty();
	}

	bool has_run_from() {
		return !run_from().empty();
	}

	// the key is the parameter's name
	std::map<std::string, Param> data;

};



#endif