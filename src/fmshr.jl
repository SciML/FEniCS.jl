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

generate_mesh(geom_object::Geometry,size::Int)=Mesh(mshr.generate_mesh(geom_object.pyobject,size))

+(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__add__](geom_object2.pyobject))
-(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__sub__](geom_object2.pyobject))
*(geom_object1::Geometry, geom_object2::Geometry) = Geometry(geom_object1.pyobject[:__mul__](geom_object2.pyobject))

export Circle,Rectangle,Ellipse,Box,Cone,Sphere,generate_mesh
