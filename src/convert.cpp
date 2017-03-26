#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>
#include <typeinfo>
#include <typeindex>
#include <cassert>
#include <unordered_map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <matrix_serialization.hpp>
#include <logging.hpp>

namespace po = boost::program_options;



enum class StreamType {Text, Bin};


template<class _OStream>
auto operator<< (_OStream & out, StreamType const& type) -> _OStream&
{
	switch(type) {
	case StreamType::Text: out << "TEXT"; break;
	case StreamType::Bin:  out << "BIN";  break;
	default: throw std::invalid_argument{"Unknown StreamType."};
	} // end switch

	return out;
}


template<class _IStream>
auto operator>> (_IStream & in, StreamType & type) -> _IStream&
{
	std::string s;
	in >> s;
	boost::to_upper(s);

	if(s == "TEXT") {
		type = StreamType::Text;
	} 
	else if(s == "BIN") {
		type = StreamType::Bin;
	} 
	else {
		throw std::invalid_argument{ "Could not convert `" 
		                           + s + "` to StreamType." };
	}

	return in;
}


using text_t = std::integral_constant<StreamType, StreamType::Text>;
using bin_t  = std::integral_constant<StreamType, StreamType::Bin>;



auto init_options() -> po::options_description
{
	po::options_description desc;

	desc.add_options()
		( "help", "Produce the help message." )
		( "type", po::value<std::string>()->required()
		, "Type of an element of the matrix. It may be one of: "
		  "float, double, cfloat, cdouble." )
		( "from", po::value<StreamType>()->required()
		, "Type of the input stream. It may be either \"Text\" or \"Bin\"."
		  " \"Text\" means a commom dat-file format: columns are separated "
		  " by TABs and rows by NEWLINEs. \"Bin\" means the binary format "
		  "of the  Boost Serialization library." )
		( "to", po::value<StreamType>()->required()
		, "Type of the output stream." );
	
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

	if (vm.count("help")) {
		help(description);
		return;
	}

	po::notify(vm);

	run( element_type(vm["type"].as<std::string>()) 
	   , vm["from"].as<StreamType>()
	   , vm["to"].as<StreamType>()
	   );
}


template<class _T, class _IStream>
auto load_line(_IStream & input, std::vector<_T> & data) -> std::size_t
{
	std::string line;
	std::getline(input, line);
	std::istringstream line_stream{line};

	auto const old_size = data.size();
	std::copy( std::istream_iterator<_T>{line_stream}
	         , std::istream_iterator<_T>{}
	         , std::back_inserter(data) );
	return data.size() - old_size;
}


template<class _T, class _IStream>
auto load( _IStream & input, text_t) -> tcm::Matrix<_T>
{
	std::vector<_T> _temp_data_;

	auto width = load_line(input, _temp_data_);

	while (input) {
		auto const current_width = load_line(input, _temp_data_);
		if (current_width and width != current_width) {
			throw std::runtime_error
				{ "Data has invalid format: width of the first row (" 
			    + std::to_string(width) + ") does not match the width of the "
			    " current row (" + std::to_string(current_width) + ")." };
		}
	}

	auto const height =	_temp_data_.size() / width;
	assert(_temp_data_.size() % width == 0);

	tcm::Matrix<_T> A{height, width};
	for(std::size_t i = 0; i < height; ++i)
		for(std::size_t j = 0; j < width; ++j)
			A(i, j) = _temp_data_[j + width * i];

	return A;
}


template<class _T, class _IStream>
auto load( _IStream & input, bin_t ) -> tcm::Matrix<_T>
{
	boost::archive::binary_iarchive input_archive{input};

	tcm::Matrix<_T> A;
	input_archive >> A;

	return A;
}


template<class _T, class _IStream, class _Logger>
auto load(_IStream & input, _Logger & lg, StreamType stream_type) 
	-> tcm::Matrix<_T>
{
	// LOG(lg, info) << "[*] Loading...";
	switch (stream_type) {
		case StreamType::Text: return load<_T>(input, text_t{});
		case StreamType::Bin:  return load<_T>(input, bin_t{});
		default: throw std::invalid_argument{"Unknown StreamType."};
	}
}



template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output, text_t ) -> void
{
	output << std::scientific << std::setprecision(20) << A;
}


template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output, bin_t ) -> void
{
	boost::archive::binary_oarchive output_archive{output};
	output_archive << A;
}


template<class _T, class _OStream, class _Logger>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output
         , _Logger & lg
         , StreamType stream_type ) -> void
{
	// LOG(lg, info) << "[*] Saving...";
	switch (stream_type) {
		case StreamType::Text: save(A, output, text_t{}); break;
		case StreamType::Bin:  save(A, output, bin_t{} ); break;
		default: throw std::invalid_argument{"Unknown StreamType."};
	}
}


template<class _T, class _IStream, class _OStream>
auto convert( _IStream & input,  StreamType input_type
            , _OStream & output, StreamType output_type ) -> void
{
	boost::log::sources::severity_logger<tcm::severity_level> lg;
	save( load<_T>(input, lg, input_type)
	    , output, lg, output_type );
	// LOG(lg, info) << "[+] Done!";
}


auto run(std::type_index type, StreamType from, StreamType to) -> void
{
	using convert_function = std::function<void()>;
	static std::unordered_map<std::type_index, convert_function> const 
	dispatch = {
		{ std::type_index(typeid(float))
		, [from, to]() { convert<float>(std::cin, from, std::cout, to); } },
		{ std::type_index(typeid(double))
		, [from, to]() { convert<double>(std::cin, from, std::cout, to); } },
		{ std::type_index(typeid(std::complex<float>))
		, [from, to]() { convert<std::complex<float>>
		                 (std::cin, from, std::cout, to); } },
		{ std::type_index(typeid(std::complex<double>))
		, [from, to]() { convert<std::complex<double>>
		                 (std::cin, from, std::cout, to); } }
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
