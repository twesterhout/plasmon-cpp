#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <logging.hpp>
#include <constants.hpp>
#include <dielectric_function_v2.hpp>
#include <matrix_serialization.hpp>



namespace po = boost::program_options;

auto init_options() -> po::options_description
{
	po::options_description desc;

	desc.add_options()
		( "help", "Produce the help message." )
		( "type", po::value<std::string>()->required()
		, "Type of an element of the matrix. It may be one of: "
		  "float, double, cfloat, cdouble." );
	desc.add(tcm::init_constants_options<double>());

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

	run(element_type(vm["type"].as<std::string>()), vm);
}


template<class _T, class _IStream>
auto read_positions(_IStream & input) -> std::vector<std::array<_T, 3>>
{
	std::vector<std::array<_T, 3>> positions;

	auto i = std::istream_iterator<_T>{input};
	while(i != std::istream_iterator<_T>{}) {
		const auto x = *i++;
		const auto y = *i++;
		const auto z = *i++;
		positions.push_back({x, y, z});
	}

	return positions;
}


template<class _T, class _IStream, class _OStream>
auto make_potential( _IStream & input, _OStream & output
                   , po::variables_map const& vm ) -> void
{
	using R = tcm::utils::Base<_T>;
	auto const positions = read_positions<R>(input);
	auto const constants = tcm::load_constants<R, double, std::map<std::string, R>>(vm);
	boost::log::sources::severity_logger<tcm::severity_level> lg;
	auto const V = tcm::coulomb::make<_T>( positions
	                                     , constants
	                                     , lg );
	
	boost::archive::binary_oarchive output_archive{output};
	output_archive << V;
}


auto run( std::type_index type
        , po::variables_map const& vm ) -> void
{
	using make_function_t = std::function<void()>;
	static std::unordered_map<std::type_index, make_function_t> const 
	dispatch = {
		{ std::type_index(typeid(float))
		, [&vm](){make_potential<float>(std::cin, std::cout, vm);} },
		{ std::type_index(typeid(double))
		, [&vm](){make_potential<double>(std::cin, std::cout, vm);} },
		{ std::type_index(typeid(std::complex<float>))
		, [&vm](){make_potential<std::complex<float>>(std::cin, std::cout, vm);} },
		{ std::type_index(typeid(std::complex<double>))
		, [&vm](){make_potential<std::complex<double>>(std::cin, std::cout, vm);} }
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

