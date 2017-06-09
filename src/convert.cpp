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
using text_t = std::integral_constant<StreamType, StreamType::Text>;
using bin_t  = std::integral_constant<StreamType, StreamType::Bin>;


template<class _OStream>
auto operator<< (_OStream & out, StreamType const& type) -> _OStream&
{
	switch (type) {
		case StreamType::Text: out << "TEXT"; break;
		case StreamType::Bin:  out << "BIN";  break;
		default: throw std::invalid_argument{"Unknown StreamType."};
	} 
	return out;
}

template<class _IStream>
auto operator>> (_IStream & in, StreamType & type) -> _IStream&
{
	using namespace std::string_literals;
	static std::unordered_map<std::string, StreamType> const types =
		{ { "text"s, StreamType::Text }
		, { "bin"s,  StreamType::Bin  }
		};
		
	std::string s;
	in >> s;
	boost::to_lower(s);
	type = types.at(s);
	return in;
}





auto init_options() -> po::options_description
{
	po::options_description desc;

	desc.add_options()
		( "help", "Produce the help message." )
		( "type", po::value<std::string>()->required()
		, "Element of the matrix. It may be one of: "
		  "float, double, cfloat, cdouble." )
		( "from", po::value<StreamType>()->required()
		, "Input stream type. It may be either \"Text\" or \"Bin\". "
		  "\"Text\" means a commom dat-file format, while \"Bin\" "
		  "is the binary format used by "
		  "the Boost Serialization library." )
		( "to", po::value<StreamType>()->required()
		, "Type of the output stream (see '--from')." )
		( "column", po::value<std::size_t>()
		, "Column of the matrix to print. If not specified, " 
		  "the whole matrix is printed." );
	
	return desc;
}


auto element_type(std::string input) -> std::type_index
{
	using namespace std::string_literals;
	static std::unordered_map<std::string, std::type_index> const types = 
		{ { "float"s,   typeid(float)                }
		, { "double"s,  typeid(double)               }
		, { "cfloat"s,  typeid(std::complex<float>)  }
		, { "cdouble"s, typeid(std::complex<double>) }
		};
	
	try {
		return types.at(boost::to_lower_copy(input));
	} catch (std::out_of_range & e) {
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
	run(vm);
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
auto load(_IStream & input, text_t) -> tcm::Matrix<_T>
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
	for (std::size_t i = 0; i < height; ++i) {
		for (std::size_t j = 0; j < width; ++j) {
			A(i, j) = _temp_data_[j + width * i];
		}
	}

	return A;
}


template<class _T, class _IStream>
auto load(_IStream & input, bin_t) -> tcm::Matrix<_T>
{
	boost::archive::binary_iarchive input_archive{input};
	tcm::Matrix<_T> A;
	input_archive >> A;
	return A;
}


template<class _T, class _IStream>
auto load(_IStream & input, StreamType stream_type) -> tcm::Matrix<_T>
{
	switch (stream_type) {
		case StreamType::Text: return load<_T>(input, text_t{});
		case StreamType::Bin:  return load<_T>(input, bin_t{});
		default: throw std::invalid_argument{"Unknown StreamType."};
	}
}



template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output
         , text_t ) -> void
{
	output << std::scientific << std::setprecision(20) 
	       << A;
}

template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output
         , bin_t ) -> void
{
	boost::archive::binary_oarchive output_archive{output};
	output_archive << A;
}


template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output
         , StreamType stream_type ) -> void
{
	switch (stream_type) {
		case StreamType::Text: save(A, output, text_t{}); break;
		case StreamType::Bin:  save(A, output, bin_t{} ); break;
		default: throw std::invalid_argument{"Unknown StreamType."};
	}
}


template<class _T, class _OStream>
auto save( tcm::Matrix<_T> const& A
         , _OStream & output
         , StreamType stream_type 
         , std::size_t const column ) -> void
{
	tcm::Matrix<_T> Aj{A.height(), 1};
	std::copy( A.cbegin_column(column), A.cend_column(column)
	         , Aj.begin_column(0) );
	save(Aj, output, stream_type);
}


template<class _T, class _IStream, class _OStream>
auto convert( _IStream & input,  StreamType input_type
            , _OStream & output, StreamType output_type ) -> void
{
	save( load<_T>(input, input_type)
	    , output, output_type );
}

template<class _T, class _IStream, class _OStream>
auto convert( _IStream & input,  StreamType input_type
            , _OStream & output, StreamType output_type
            , std::size_t const column ) -> void
{
	save( load<_T>(input, input_type)
	    , output, output_type, column );
}


template <class _T>
auto run(po::variables_map const& vm) -> void
{
	auto const from = vm["from"].as<StreamType>();
	auto const to   = vm["to"  ].as<StreamType>();

	if (vm.count("column")) {
		auto const column = vm["column"].as<std::size_t>();
		convert<_T>(std::cin, from, std::cout, to, column);	
	} else {
		convert<_T>(std::cin, from, std::cout, to);	
	}
}


auto dispatch(po::variables_map const& vm) -> void
{
	std::unordered_map< std::type_index, 
		void (*)(po::variables_map const&)> const callbacks = 
			{ { typeid(float), &run<float> }
			, { typeid(double), &run<double> }
			, { typeid(std::complex<float>), &run<std::complex<float>> }
			, { typeid(std::complex<double>), &run<std::complex<double>> }
			};
	callbacks.at(element_type(vm["type"].as<std::string>()))(vm);
}



int main(int argc, char** argv)
{
	process_command_line
		( argc, argv
		, [](auto desc) { std::cout << desc << '\n'; }
		, &dispatch
		);

	return 0;
}
