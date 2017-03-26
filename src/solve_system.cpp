#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <typeinfo>
#include <typeindex>
#include <cassert>
#include <unordered_map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <logging.hpp>
#include <matrix_serialization.hpp>
#include <lapack.hpp>



namespace po = boost::program_options;



auto init_options() -> po::options_description
{
	po::options_description desc;

	desc.add_options()
		( "help", "Produce the help message." )
		( "type", po::value<std::string>()->required()
		, "Type of an element of the hamiltonian matrix. It may be one of: "
		  "float, double, cfloat, cdouble." )
		( "energies", po::value<std::string>()->required()
		, "Name of the file where to save the eigenenergies of the system." )
		( "states", po::value<std::string>()->required()
		, "Name of the file where to save the eigenestates of the system." );
	
	return desc;
}


auto element_type(std::string input) -> std::type_index
{
	using namespace std::string_literals;
	static std::unordered_map<std::string, std::type_index> const types = 
		{ { "float"s,   std::type_index(typeid(float))                }
		, { "double"s,  std::type_index(typeid(double))               }
		, { "cfloat"s,  std::type_index(typeid(std::complex<float>))  }
		, { "cdouble"s, std::type_index(typeid(std::complex<double>)) }
		};
	
	boost::to_lower(input);
	try {
		return types.at(input);
	} catch(std::out_of_range & e) {
		std::cerr << "Invalid element type `" + input + "`!\n";
		throw;
	}
}



template<class _Help, class _Run>
auto process_command_line( int argc, char** argv
                         , _Help&& help
						 , _Run&& run ) -> void
{
	auto const description = init_options();
	po::variables_map vm;

	po::store(po::command_line_parser(argc, argv).options(description).run()
	         , vm );

	if(vm.count("help")) {
		help(description);
		return;
	}

	po::notify(vm);

	run( element_type(vm["type"].as<std::string>()) 
	   , vm["energies"].as<std::string>()
	   , vm["states"].as<std::string>()
	   );
}


template<class _T, class _IStream, class _OStream1, class _OStream2>
auto solve( _IStream & input
          , _OStream1 & energies_output
		  , _OStream2 & states_output ) -> void
{
	boost::log::sources::severity_logger<tcm::severity_level> lg;

	LOG(lg, info) << "Reading Hamiltonian...";
	tcm::Matrix<_T> H;
	boost::archive::binary_iarchive input_archive{input};
	input_archive >> H;

	tcm::Matrix<tcm::utils::Base<_T>>   E{H.height(), 1};
	tcm::Matrix<_T>                   Psi{H.height(), H.height()};

	LOG(lg, info) << "Diagonalizing...";
	tcm::lapack::heevr(H, E, Psi);

	LOG(lg, info) << "Saving results...";
	boost::archive::binary_oarchive energies_archive{energies_output};
	energies_archive << E;
	boost::archive::binary_oarchive states_archive{states_output};
	states_archive << Psi;

	LOG(lg, info) << "Done!";
}


auto run( std::type_index type
        , std::string const& energies_filename
		, std::string const& states_filename ) -> void
{
	std::ofstream energies_file{energies_filename};
	if(not energies_file) {
		throw std::runtime_error{ "Could not open `" + energies_filename
		                        + "` for writing." };
	}
	std::ofstream states_file{states_filename};
	if(not states_file) {
		throw std::runtime_error{ "Could not open `" + states_filename
		                        + "` for writing." };
	}

	using solve_function_t = std::function<void()>;
	static std::unordered_map<std::type_index, solve_function_t> const 
	dispatch = {
		{ std::type_index(typeid(float))
		, [&energies_file, &states_file]() 
		  {solve<float>(std::cin, energies_file, states_file); } },
		{ std::type_index(typeid(double))
		, [&energies_file, &states_file]()
		  {solve<double>(std::cin, energies_file, states_file); } },
		{ std::type_index(typeid(std::complex<float>))
		, [&energies_file, &states_file]()
		  {solve<std::complex<float>>(std::cin, energies_file, states_file); } },
		{ std::type_index(typeid(std::complex<double>))
		, [&energies_file, &states_file]() 
		  {solve<std::complex<double>>(std::cin, energies_file, states_file);} }
	};

	tcm::setup_console_logging();
	dispatch.at(type)();
}



int main(int argc, char** argv)
{
	process_command_line
		( argc, argv
		, [](auto desc) { std::cout << desc << '\n'; }
		, &run 
		);
	return 0;
}

