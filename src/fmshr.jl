@pyimport mshr
@fenicsclass Geometry

#functions necessary for creating meshes from geometrical objects.
#2D objects below
Circle(centre,radius) = Geometry(mshr.Circle(centre,radius))
Rectangle(corner1,corner2)=Geometry(mshr.Rectangle(corner1,corner2))
Ellipse(centre,horizontal_semi_axis,vertical_semi_axis,fragments)=Geometry(mshr.Ellipse(centre,horizontal_semi_axis,vertical_semi_axis,fragments))

#3d objects below
Box(corner1,corner2) = Geometry(mshr.Box(corner1,corner2))
Cone(top,bottom,bottom_radius,slices::Int)=Geometry(mshr.Cone(top,bottom,bottom_radius,slices))
Sphere(centre,radius) = Geometry(mshr.Sphere(centre,radius))

#we generate the mesh of some geometrical objects
generate_mesh(geom_object::Geometry,size::Int)=Mesh(mshr.generate_mesh(geom_object.pyobject,size))
export Circle,Rectangle,Ellipse,Box,Cone,Sphere,generate_mesh


#operator overloading for Geometry types so we can create composite shapes
+(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__add__](geom_object2.pyobject))
-(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__sub__](geom_object2.pyobject))
*(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__mul__](geom_object2.pyobject))

set_subdomain(object::Geometry,degree,domain2) = fenicspycall(object, :set_subdomain,degree,domain2.pyobject)
"""
Extrude2D

Extrudes a 2D geometry to 3D
"""
Extrude2D(object::Geometry,thickness) = Geometry(mshr.Extrude2D(object.pyobject,thickness))
export  set_subdomain, Extrude2D

CSGUnion(geom_object1::Geometry, geom_object2::Geometry) = Geometry(mshr.CSGUnion(geom_object1.pyobject,geom_object2.pyobject))
CSGIntersection(geom_object1::Geometry, geom_object2::Geometry) = Geometry(mshr.CSGIntersection(geom_object1.pyobject,geom_object2.pyobject))
CSGDifference(geom_object1::Geometry, geom_object2::Geometry) = Geometry(mshr.CSGDifference(geom_object1.pyobject,geom_object2.pyobject))

CSGScaling(geom_object1::Geometry, scale) = Geometry(mshr.CSGScaling(geom_object1.pyobject,scale))

CSGRotation(geom_object1::Geometry, angle) = Geometry(mshr.CSGRotation(geom_object1.pyobject,angle))
CSGRotation(geom_object1::Geometry, point, angle) = Geometry(mshr.CSGRotation(geom_object1.pyobject,point, angle))
CSGRotation(geom_object1::Geometry, point_axis, point_center, angle) = Geometry(mshr.CSGRotation(geom_object1.pyobject, point_axis, point_center, angle))


CSGTranslation(geom_object1::Geometry, translation_point) = Geometry(mshr.CSGTranslation(geom_object1.pyobject, translation_point))

export CSGUnion, CSGIntersection, CSGDifference, CSGScaling, CSGRotation, CSGTranslation
