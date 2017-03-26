#ifndef TCM_IO_HPP
#define TCM_IO_HPP

#include <fstream>
#include <vector>
#include <algorithm>

#include <matrix.hpp>
#include <logging.hpp>
#include <benchmark.hpp>



namespace tcm {


namespace io {

template<class T, class Logger>
auto hamiltonian_from_text( std::string const& filename
                          , Logger & lg ) -> Matrix<T>
{
	MEASURE;
	LOG(lg, debug) << "Reading hamiltonian from text file...";

	std::ifstream in_stream{filename};
	if(!in_stream) {
		LOG(lg, error) << "Could not open `" << filename << "` for reading";
		throw std::runtime_error{"Could not open file."};
	}

	std::vector<T> data;
	std::copy( std::istream_iterator<T>{in_stream}, std::istream_iterator<T>{}
	         , std::back_inserter(data) );
	auto const N = static_cast<std::size_t>(std::round(std::sqrt(data.size())));

	Matrix<T> H{N, N};
	std::copy(std::begin(data), std::end(data), H.data());
	
	LOG(lg, debug) << "Successfully read the hamiltonian...";
	return H;
}


template<class T, class Logger>
auto hamiltonian_from_bin( std::string const& filename
                         , Logger & lg ) -> Matrix<T>
{
	MEASURE;
	LOG(lg, info) << "Reading hamiltonian from binary file...";

	std::ifstream in_stream{filename};
	if(!in_stream) {
		LOG(lg, error) << "Could not open `" << filename << "` for reading";
		throw std::runtime_error{"Could not open file."};
	}

	Matrix<T> H;
	boost::archive::binary_iarchive in_archive{in_stream};
	in_archive >> H;
	
	LOG(lg, debug) << "Successfully read the hamiltonian...";
	return H;
}



//! \brief Reads cite positions from an input stream.

//! This function assumes that the input consists of three
//! columns representing the x, y, z coordinates of the cites.
template<class T, class Logger>
auto positions_from_text( std::string const& filename
                        , Logger & lg ) -> std::vector<std::array<T, 3>>
{
	MEASURE;
	LOG(lg, debug) << "Reading atomic positions from text file...";

	std::ifstream in_stream{filename};
	if(!in_stream) {
		LOG(lg, error) << "Could not open `" << filename << "` for reading";
		throw std::runtime_error{"Could not open file."};
	}

	std::vector< std::array<T, 3> > positions;
	auto i = std::istream_iterator<T>{in_stream};
	while(i != std::istream_iterator<T>{}) {
		const auto x = *i++;
		const auto y = *i++;
		const auto z = *i++;
		positions.push_back({x, y, z});
	}

	LOG(lg, debug) << "Successfully read atomic positions...";
	return positions;
}



} // namespace io


} // namespace tcm


#endif // TCM_IO_HPP_
