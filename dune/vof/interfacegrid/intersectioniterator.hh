#ifndef DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH
#define DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH

#include <dune/vof/interfacegrid/intersection.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIntersectionIterator
    // --------------------------

    template< class Grid, class HostIntersectionIterator >
    class InterfaceGridIntersectionIterator
    {
    protected:
      typedef InterfaceGridIntersectionIterator< Grid, HostIntersectionIterator > This;

      typedef typename remove_const< Grid >::type::Traits Traits;

      static const bool isLeafIntersection =
        is_same< HostIntersectionIterator,
                 typename Grid::HostGrid::Traits::LeafIntersectionIterator > :: value ;
    public:
      typedef typename conditional< isLeafIntersection,
                                    typename Traits :: LeafIntersection,
                                    typename Traits :: LevelIntersection > :: type  Intersection ;
      typedef typename Intersection :: Implementation IntersectionImpl ;

      typedef typename Traits :: ExtraData ExtraData;

      InterfaceGridIntersectionIterator ()
       : hostIterator_(),
         data_()
      {}

      InterfaceGridIntersectionIterator ( ExtraData data, const HostIntersectionIterator &hostIterator )
       : hostIterator_( hostIterator ),
         data_( data )
      {}

      InterfaceGridIntersectionIterator ( const This &other )
       : hostIterator_( other.hostIterator_ ),
         data_( other.data() )
      {}

      bool equals ( const This &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }

      void increment ()
      {
        ++hostIterator_;
      }

      Intersection dereference () const
      {
        return IntersectionImpl( data(), *hostIterator_ );
      }

      ExtraData data() const { return data_; }

    protected:
      HostIntersectionIterator hostIterator_;
      ExtraData data_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH
