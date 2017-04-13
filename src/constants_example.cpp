#include <iostream>
#include <constants.hpp>

namespace po = boost::program_options;

auto init_options() -> po::options_description
{
	po::options_description description{"General options"};
	// Add some general options
	description.add_options()
		( "help", "Produce the help message." )
		( "in.constants.hello-world"
		, boost::program_options::value<double>()->default_value(123456789.)
		, "Dummy value" );
	// Now add options to specify physical constants.
	description.add(tcm::init_constants_options<double>());
	return description;
}


int main(int argc, char** argv)
{
	auto const desc = init_options();
	po::variables_map vm;
	po::store( po::command_line_parser(argc, argv)
	             .options(desc)
	             .run()
	         , vm );

	if (vm.count("help")) {
		std::cout << desc << '\n';
		return 0;
	}

	po::notify(vm);

	// Important that in.constants.hello-world is also loaded!
	// Not only the options from tcm::default_constants().
	auto const cs = tcm::load_constants<double>(vm);

	for_each( std::begin(cs), std::end(cs)
	        , [](auto const& i) { std::cout << i.first 
	                                        << " = " 
	                                        << i.second << '\n'; } );
	return 0;
}
