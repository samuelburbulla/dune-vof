#ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
#define DUNE_VOF_INTERFACEGRID_INTERSECTION_HH

#include <dune/vof/interfacegrid/declaration.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIntersection
    // ------------------

    template< class Grid, class HostIntersection >
    class InterfaceGridIntersection
    {
    protected:
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits :: ExtraData ExtraData ;

    public:
      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity         Entity;
      typedef typename Traits::template Codim< 0 >::EntityImpl     EntityImpl;
      typedef typename Traits::template Codim< 1 >::Geometry       Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry  LocalGeometry;

    public:
      InterfaceGridIntersection ()
      : hostIntersection_(),
        data_()
      {}

      explicit InterfaceGridIntersection ( ExtraData data )
      : hostIntersection_(),
        data_( data )
      {}

      InterfaceGridIntersection ( ExtraData data, const HostIntersection &hostIntersection )
      : hostIntersection_( hostIntersection ),
        data_( data )
      {}

      // template the other intersection here since it could be leaf or level
      // intersection and we don't want to specify this here
      template < class IntersectionImpl >
      bool equals ( const IntersectionImpl& other ) const
      {
        return hostIntersection() == other.hostIntersection();
      }

      Entity inside () const
      {
        return Entity( EntityImpl( data(), hostIntersection().inside() ) );
      }

      Entity outside () const
      {
        return Entity( EntityImpl( data(), hostIntersection().outside() ) );
      }

      bool boundary () const { return hostIntersection().boundary(); }

      bool conforming () const { return hostIntersection().conforming(); }

      bool neighbor () const { return hostIntersection().neighbor(); }

      int boundaryId () const { return hostIntersection().boundaryId(); }

      size_t boundarySegmentIndex () const
      {
        return hostIntersection().boundarySegmentIndex();
      }

      LocalGeometry geometryInInside () const
      {
        return LocalGeometry( hostIntersection().geometryInInside() );
      }

      LocalGeometry geometryInOutside () const
      {
        return LocalGeometry( hostIntersection().geometryInOutside() );
      }

      Geometry geometry () const
      {
        return Geometry( hostIntersection().geometry() );
      }

      GeometryType type () const { return hostIntersection().type(); }

      int indexInInside () const { return hostIntersection().indexInInside(); }
      int indexInOutside () const { return hostIntersection().indexInOutside(); }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return hostIntersection().integrationOuterNormal( local );
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return hostIntersection().outerNormal( local );
      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return hostIntersection().unitOuterNormal( local );
      }

      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        return hostIntersection().centerUnitOuterNormal();
      }

      const HostIntersection &hostIntersection () const
      {
        return hostIntersection_;
      }

      ExtraData data() const { return data_; }

    protected:
      HostIntersection hostIntersection_;
      ExtraData  data_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
