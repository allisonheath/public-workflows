/*
    libmaus
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#if ! defined(LIBMAUS_BAMBAM_BAMMERGETEMPLATE_HPP)
#define LIBMAUS_BAMBAM_BAMMERGETEMPLATE_HPP

#include <libmaus/bambam/BamCatHeader.hpp>

namespace libmaus
{
	namespace bambam
	{
		/**
		 * class for merging BAM input
		 **/
		template<typename _heap_comparator_type, typename _sort_check_type>
		struct BamMergeTemplate : public BamAlignmentDecoder, public BamAlignmentDecoderWrapper
		{
			private:
			typedef _heap_comparator_type heap_comparator_type;
			typedef _sort_check_type sort_check_type;
			typedef BamMergeTemplate<heap_comparator_type,sort_check_type> this_type;
			typedef typename libmaus::util::unique_ptr<this_type>::type unique_ptr_type;
			
			std::vector<libmaus::bambam::BamAlignmentDecoderInfo> infos;
			std::vector<std::string> const filenames;
			libmaus::bambam::BamCatHeader header;
			libmaus::autoarray::AutoArray<libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type> wrappers;
			libmaus::autoarray::AutoArray<libmaus::bambam::BamAlignmentDecoder *> decoders;
			libmaus::autoarray::AutoArray<libmaus::bambam::BamAlignment *> algns;
			heap_comparator_type comp;
			std::priority_queue < uint64_t, std::vector< uint64_t >, heap_comparator_type > Q;
			
			void tryLoad(uint64_t const id)
			{
				if ( decoders[id]->readAlignment() )
				{
					header.updateAlignment(id,*algns[id]);
					Q.push(id);
				}
			}
			
			void init()
			{
				if ( ! sort_check_type::issorted(header) )
				{
					libmaus::exception::LibMausException se;
					se.getStream() << "BamMergeTemplate::BamMergeTemplate(): cannot merge, not all files are marked as sorted." << std::endl;
					se.finish();
					throw se;
				}
				
				if ( ! sort_check_type::istopological(header) )
				{
					libmaus::exception::LibMausException se;
					se.getStream() << "BamMergeTemplate::BamMergeTemplate(): cannot merge, order of reference sequences is inconsistent between files (no topological sorting of order graph possible)." << std::endl;
					se.finish();
					throw se;				
				}
				
				header.bamheader->changeSortOrder(sort_check_type::getSortOrder());
			
				for ( uint64_t i = 0; i < infos.size(); ++i )
				{
					libmaus::bambam::BamAlignmentDecoderWrapper::unique_ptr_type tptr ( libmaus::bambam::BamAlignmentDecoderFactory::construct(infos[i]) );
					wrappers[i] = UNIQUE_PTR_MOVE(tptr);
					decoders[i] = &(wrappers[i]->getDecoder());
					algns[i] = &(decoders[i]->getAlignment());
					tryLoad(i);
				}
			}
			
			public:
			BamMergeTemplate(std::vector<std::string> const & rfilenames, bool const putrank = false) 
			: BamAlignmentDecoder(putrank), 
			  infos(libmaus::bambam::BamAlignmentDecoderInfo::filenameToInfo(rfilenames)),
			  header(infos),
			  wrappers(infos.size()),
			  decoders(infos.size()), 
			  algns(infos.size()), 
			  comp(algns.begin()), 
			  Q(comp)
			{
				init();
			}

			BamMergeTemplate(std::vector<libmaus::bambam::BamAlignmentDecoderInfo> const & rinfos, bool const putrank = false)
			: BamAlignmentDecoder(putrank), 
			  infos(rinfos),
			  header(infos),
			  wrappers(infos.size()),
			  decoders(infos.size()), 
			  algns(infos.size()), 
			  comp(algns.begin()), 
			  Q(comp)
			{
				init();
			}

			/**
			 * input method
			 *
			 * @return bool if alignment input was successfull and a new alignment was stored
			 **/
			virtual bool readAlignmentInternal(bool const delayPutRank = false)
			{
				if ( expect_true(!Q.empty()) )
				{
					uint64_t const id = Q.top();
					Q.pop();
					
					alignment.swap(*algns[id]);

					if ( ! delayPutRank )
						putRank();

					tryLoad(id);
					
					return true;
				}
				else
				{
					return false;
				}
			}
			
			virtual libmaus::bambam::BamHeader const & getHeader() const
			{
				return *(header.bamheader);
			}
			
			virtual libmaus::bambam::BamAlignmentDecoder & getDecoder()
			{
				return *this;
			}
		};
	}
}
#endif
