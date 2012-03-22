<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex12",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 12 - The <code>MeshData</code> class</h1>

<br><br>The previous examples covered the certainly involved
aspects of simulating multiple equation systems, and prior 
to this even using adaptive refinement.  I.e.: cool stuff.

<br><br>This example now reduces brain effort significantly,
and presents some supplements concerning data I/O,
especially how to handle the <code> MeshData </code> object.
The <code> MeshData </code> class may be used for handling input
data for a simulation, like actual material properties, 
(not indicators, like in the <code> BoundaryInfo </code> class),
or mode shapes that serve as input for acoustic radiation
simulations.  The use of the <code> MeshData </code> during simulation
is straightforward:  the 
<ul>
<li>
<code>Number MeshData::operator()(const Node*, int)</code>,
get the i-th floating-point value associated with the node
</li>
<li>
<code> bool MeshData::has_data(const Node*)</code>,
verify whether a certain node has data associated
</li>
<li>
<code>std::vector<Number>& MeshData::get_data (const Node* node)</code>
to get read-access to the data associated with this node (better
make sure first that this node <i>has</i> data, see <code> has_data() )
</li>
<li>
iterator for nodal data <code>MeshData::const_node_data_iterator</code>
to directly iterate over the set of nodes that hold data, 
instead of asking through <code>has_data()</code> each time again.
</li>
</ul>
(and corresponding methods for <code> const Elem*</code>) provide access to
the floating-point data associated with nodes/elements.
This example does <i>not</i> show how to use these aforementioned
methods, this is straightforward.  Simply check if the current 
<code>Elem*</code> has a node with data. If so, add this contribution to the 
RHS, or whatever. Or ask the <code> MeshData </code> for each <code>Elem*</code> for its 
material properties...

<br><br>Here, only the actual file I/O necessary to handle such 
nodal/element data is presented.  Lean back, relax...


<br><br>C++ include files that we need.
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Functions to initialize the library
and provide some further features (e.g. 
our own pi)
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh and 
<code> MeshData </code> functionality.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        #include "mesh_tools.h"
        #include "mesh_data.h"
        #include "unv_io.h"
        #include "gmv_io.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric vertex associated with a Mesh
</div>

<div class ="fragment">
<pre>
        #include "node.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Function prototype for creating artificial nodal data
that can be inserted into a <code>MeshData</code> object.
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh& mesh,
                                     std::map&lt;const Node*, std::vector&lt;Number&gt; &gt;& art_data);
        
</pre>
</div>
<div class = "comment">
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize the library.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
          if (libMesh::n_processors() &gt; 1)
            {
              if (libMesh::processor_id() == 0)
                {
                  std::cerr &lt;&lt; "ERROR: Skipping example 12. " &lt;&lt; std::endl;
                  std::cerr &lt;&lt; "MeshData objects currently only work in serial." &lt;&lt; std::endl;
                }
              return 0;
            }
        
</pre>
</div>
<div class = "comment">
Check for proper usage. The program is designed to be run
as follows:
<pre>
$ ./ex12 -d 3 in_mesh.unv
</pre>
where in_mesh.unv should be a Universal file.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 4)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -d &lt;dim&gt; in_mesh.unv"
                          &lt;&lt; std::endl;
              
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
Get the dimensionality of the mesh from argv[2]
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = std::atoi(argv[2]);
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D/3D support");
          
</pre>
</div>
<div class = "comment">
The filename of the mesh
</div>

<div class ="fragment">
<pre>
          const std::string mesh_file = argv[3];
        
</pre>
</div>
<div class = "comment">
The following example makes currently sense
only with a Universal file
</div>

<div class ="fragment">
<pre>
          if (mesh_file.rfind(".unv") &gt;= mesh_file.size())
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "ERROR:  This example works only properly with a Universal mesh file!"
                        &lt;&lt; std::endl;
              
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
Some notes on <code>MeshData</code>:
<ul>
<li>
The <code>MeshData</code> is <i>not</i> a mesh!  Consult the <code>Mesh</code>,
<code>MeshBase</code>, and <code>BoundaryMesh</code> classes for details on
handling meshes!
</li>
<li>
The <code>MeshData</code> is an object handling arbitrary floating-point 
data associated with mesh entities (nodes, elements), currently
most features only available for nodal data.  However,
extending to element-data (does <i>anybody</i> need this?)
is straightforward.
</li>
<li>
Currently std::vector<Number> is the data (associated
with nodes/elements) that can be handled by <code>MeshData</code>,
</li>
<li>
In order to provide <i>full</i> functionality, the <code>MeshData</code>
<i>has</i> to be activated prior to reading in a mesh.  However,
there is a so-called compatibility mode available when you 
(intentionally) forgot to "turn on" the <code>MeshData</code>.
</li>
<li>
It is possible to provide user-defined nodal data that
may subsequently be used or written to file.
</li>
<li>
Translate the nodal-/element-associated data to formats
suitable for visual inspection, e.g. GMV files.
</li>
<li>
Two I/O file formats are currently supported: the Universal 
file format with dataset 2414, and an extension of 
the libMesh-own xda/xdr format, named xtr, xta.
Some details on this:
</li>
<li>
The xtr/xta format is simply an extension of the
xdr/xda format to read/write also nodal/element-associated
floating-point data.  The xtr interface of the <code>MeshData</code>
creates stand-alone files that may be binary or ASCII.
You cannot append files created by the <code>MeshData</code> I/O methods
to a mesh file (xdr/xda).
The xtr files may be seen as a "backdrop", especially when
binary files are preferred.  Note that unv files are <i>always</i>
ASCII and may become pretty big!
</li>
<li>
The unv format is an extremely versatile text-based file format
for arbitrary data.  Its functionality is <i>large</i>, and <code>libMesh</code>
supports only a small share: namely the I/O for nodes, elements and
arbitrary node-/element-associated data, stored in the 
so-called datasets "2411", "2412", and "2414", respectively.
Here, only the last dataset is covered.  The two first datasets are
implemented in the Universal support for I/O of meshes.  A single
unv file may hold <i>multiple</i> datasets of type "2414".  To
distinguish data, each dataset has a header.  The <code>libMesh</code> pendant
to this header is the <code>MeshDataUnvHeader</code> class.  When you
read/write unv files using the <code>MeshData</code>, you <i>always</i>
automatically also read/write such headers.
</li>
</ul>
Enough babble for now.  Examples. 
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Create a mesh with the requested dimension.
Below we require a 3-dim mesh, therefore assert
it.
</div>

<div class ="fragment">
<pre>
            libmesh_assert (dim == 3);
            Mesh mesh;
            MeshData mesh_data(mesh);
          
</pre>
</div>
<div class = "comment">
Activate the <code>MeshData</code> of the mesh, so that
we can remember the node and element ids used
in the file.  When we do not activate the <code>MeshData</code>,
then there is no chance to import node- or element-
associated data.
</div>

<div class ="fragment">
<pre>
            mesh_data.activate();
        
</pre>
</div>
<div class = "comment">
Now we can safely read the input mesh.  Note that
this should be an .xda/.xdr or .unv file only, otherwise
we cannot load/save associated data.
</div>

<div class ="fragment">
<pre>
            mesh.read (mesh_file, &mesh_data);
            
</pre>
</div>
<div class = "comment">
Print information about the mesh and the data
to the screen.  Obviously, there is no data
(apart from node & element ids) in it, yet.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "Finished reading the mesh.  MeshData is active but empty:" &lt;&lt; std::endl
                      &lt;&lt; "---------------------------------------------------------" &lt;&lt; std::endl;
            
            mesh.print_info();
            mesh_data.print_info();
          
</pre>
</div>
<div class = "comment">
Create some artificial node-associated data and store
it in the mesh's <code>MeshData</code>.  Use a <code>std::map</code> for this. 
Let the function <code>create_artificial_data()</code> do the work.
</div>

<div class ="fragment">
<pre>
            {
              std::map&lt;const Node*, std::vector&lt;Number&gt; &gt; artificial_data;
              
              create_artificial_data (mesh, artificial_data);
              
</pre>
</div>
<div class = "comment">
Before we let the map go out of scope, insert it into
the <code>MeshData</code>.  Note that by default the element-associated
data containers are closed, so that the <code>MeshData</code> is
ready for use.
</div>

<div class ="fragment">
<pre>
              mesh_data.insert_node_data(artificial_data);
        
</pre>
</div>
<div class = "comment">
Let <code>artificial_data()</code> go out of scope
</div>

<div class ="fragment">
<pre>
            }
            
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "After inserting artificial data into the MeshData:" &lt;&lt; std::endl
                      &lt;&lt; "--------------------------------------------------" &lt;&lt; std::endl;
            
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Write the artificial data into a universal
file.  Use a default header for this.
</div>

<div class ="fragment">
<pre>
            std::string first_out_data="data_first_out_with_default_header.unv";
            std::cout &lt;&lt; "Writing MeshData to: " &lt;&lt; first_out_data &lt;&lt; std::endl;
            mesh_data.write(first_out_data);
        
</pre>
</div>
<div class = "comment">
Alternatively, create your own header.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "Attach our own MeshDataUnvHeader to the MeshData:" &lt;&lt; std::endl
                      &lt;&lt; "-------------------------------------------------" &lt;&lt; std::endl
                      &lt;&lt; " (note the warning: the number of values per node in" &lt;&lt; std::endl
                      &lt;&lt; "  my_header is not correct)" &lt;&lt; std::endl &lt;&lt; std::endl;
            MeshDataUnvHeader my_header;
        
</pre>
</div>
<div class = "comment">
Specify an int for this dataset.
This is particularly helpful
when there are multiple datasets in
a file and you want to find a specific one.
For this, use MeshDataUnvHeader::which_dataset(int),
which is not covered here.
</div>

<div class ="fragment">
<pre>
            my_header.dataset_label = 3;
        
</pre>
</div>
<div class = "comment">
Specify some text that helps the <i>user</i> 
identify the data.  This text is <i>not</i> used for
finding a specific dataset.  Leave default for
the remaining 2 lines.
</div>

<div class ="fragment">
<pre>
            my_header.id_lines_1_to_5[0] = "Artificial data";
            my_header.id_lines_1_to_5[1] = "sin curve in z-direction";
            my_header.id_lines_1_to_5[2] = "line in x-direction";
            
</pre>
</div>
<div class = "comment">
Also some float data, not associated with nodes,
can be stored in the header.
</div>

<div class ="fragment">
<pre>
            my_header.record_12[0] = libMesh::pi;
            
</pre>
</div>
<div class = "comment">
Now attach this header to the <code>MeshData</code>, and write
the same file again, but with the personalized header.
</div>

<div class ="fragment">
<pre>
            mesh_data.set_unv_header(&my_header);
        
</pre>
</div>
<div class = "comment">
Write again to file.
</div>

<div class ="fragment">
<pre>
            std::string second_out_data="data_second_with_header_out.unv";
            std::cout &lt;&lt; "Writing MeshData to: " &lt;&lt; second_out_data &lt;&lt; std::endl;
            mesh_data.write(second_out_data);
            
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "Before clearing the MeshData:" &lt;&lt; std::endl
                      &lt;&lt; "-----------------------------" &lt;&lt; std::endl;
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Now clear only the data associated with nodes/elements,
but keep the node/element ids used in the mesh.  To
clear also these ids, use MeshData::slim() (not used here).
</div>

<div class ="fragment">
<pre>
            mesh_data.clear();
        
</pre>
</div>
<div class = "comment">
Print information about the data to the screen.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "After clearing the MeshData:" &lt;&lt; std::endl
                      &lt;&lt; "----------------------------" &lt;&lt; std::endl;
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Now the <code>MeshData</code> is open again to read data from
file.  Read the file that we created first.
</div>

<div class ="fragment">
<pre>
            mesh_data.read(first_out_data);
            
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "After re-reading the first file:" &lt;&lt; std::endl
                      &lt;&lt; "--------------------------------" &lt;&lt; std::endl;
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Let the mesh, the unv_header etc go out of scope, and
do another example.
</div>

<div class ="fragment">
<pre>
          }
        
        
            
          std::cout &lt;&lt; std::endl 
                    &lt;&lt; "----------------------------------------------" &lt;&lt; std::endl
                    &lt;&lt; "---------- next example with MeshData --------" &lt;&lt; std::endl
                    &lt;&lt; "----------------------------------------------" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Create a new mesh, read it again -- but this time
with de-activated <code>MeshData</code>.  Then we are
able to use the compatibility mode not only for
handling <code>MeshData</code> dat, but also to <i>write</i> 
a mesh in unv format.
The libMesh-internal node and element ids are used.
</div>

<div class ="fragment">
<pre>
          {
            Mesh mesh(dim);
            MeshData mesh_data(mesh);
            
</pre>
</div>
<div class = "comment">
Read the input mesh, but with deactivated <code>MeshData</code>.
</div>

<div class ="fragment">
<pre>
            mesh.read(mesh_file, &mesh_data);
            
</pre>
</div>
<div class = "comment">
Print information about the mesh and the data
to the screen.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "De-activated MeshData:" &lt;&lt; std::endl
                      &lt;&lt; "----------------------" &lt;&lt; std::endl;
            mesh.print_info();
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Write the <i>mesh</i> (not the MeshData!) as .unv file.
In general, the <code>MeshBase</code> interface for .unv I/O
needs an active <code>MeshData</code>.  However, use compatibility
mode to at least write a .unv file, but with the 
ids from libMesh.
</div>

<div class ="fragment">
<pre>
            const std::string out_mesh = "mesh_with_libmesh_ids.unv";
            std::cout &lt;&lt; "Writing _Mesh_ to: " &lt;&lt; out_mesh &lt;&lt; std::endl
                      &lt;&lt; "Try 'diff " &lt;&lt; out_mesh &lt;&lt; " " &lt;&lt; mesh_file &lt;&lt; "'" &lt;&lt; std::endl
                      &lt;&lt; "to see the differences in node numbers." &lt;&lt; std::endl
                      &lt;&lt; "---------------------------------------" &lt;&lt; std::endl
                      &lt;&lt; std::endl;
            mesh.write(out_mesh, &mesh_data);
        
</pre>
</div>
<div class = "comment">
Again create some artificial node-associated data,
as before.
</div>

<div class ="fragment">
<pre>
            {
              std::map&lt;const Node*, std::vector&lt;Number&gt; &gt; artificial_data;
              create_artificial_data (mesh, artificial_data);
              mesh_data.insert_node_data(artificial_data);
            }
        
</pre>
</div>
<div class = "comment">
Note that even with (only) compatibility mode MeshData
data can be written.  But again, the ids from libMesh
are used.  Consult the warning messages issued in
DEBUG mode.  And the user <i>has</i> to specify that the
<code>MeshData</code> should change to compatibility mode.
</div>

<div class ="fragment">
<pre>
            mesh_data.enable_compatibility_mode();
        
</pre>
</div>
<div class = "comment">
Now that compatibility mode is used, data can be written.
_Without_ explicitly enabling compatibility mode, we
would get an error message and (with the following write()
statement).
</div>

<div class ="fragment">
<pre>
            std::string mesh_data_file = "data_third_with_libmesh_ids_out.unv";
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "Writing MeshData to: " &lt;&lt; mesh_data_file &lt;&lt; std::endl
                      &lt;&lt; "----------------------------------------------------------" 
                      &lt;&lt; std::endl &lt;&lt; std::endl;
            mesh_data.write (mesh_data_file);
        
        #ifdef HAVE_ZLIB_H
        
</pre>
</div>
<div class = "comment">
As may already seen, UNV files are text-based, so they may
become really big.  When <code>./configure</code> found <code>zlib.h</code>,
then we may also <i>read</i> or <i>write</i> <code>.unv</code> files in gzip'ed 
format! -- Pretty cool, and also pretty fast, due to zlib.h.

<br><br>Note that this works also for mesh files, not only for 
meshdata files.

<br><br>In order to write a ".unv.gz" file instead of a ".unv" file,
simply provide the full name with ".gz" appended to the
write method; it will then figure out whether this file should
be gzip'ed or not.
</div>

<div class ="fragment">
<pre>
            std::string packed_mesh_data_file =  "packed_" + mesh_data_file + ".gz";
            std::cout &lt;&lt; std::endl 
                      &lt;&lt; "Writing gzip'ed MeshData to: " &lt;&lt; packed_mesh_data_file &lt;&lt; std::endl
                      &lt;&lt; "---------------------------------------------------------------------------" &lt;&lt; std::endl
                      &lt;&lt; " To verify the integrity of the packed version, type:" &lt;&lt; std::endl &lt;&lt; std::endl
                      &lt;&lt; "   gunzip " &lt;&lt; packed_mesh_data_file &lt;&lt; "; " &lt;&lt; std::endl
                      &lt;&lt; "   diff packed_" &lt;&lt; mesh_data_file &lt;&lt; " " 
                      &lt;&lt; mesh_data_file &lt;&lt; std::endl &lt;&lt; std::endl;
            
            mesh_data.write (packed_mesh_data_file);
            
        #endif
        
        
              
</pre>
</div>
<div class = "comment">
And now a last gimmick: The <code>MeshData::translate()</code>
conveniently converts the nodal- or element-associated
data (currently only nodal) to vectors that may be used 
for writing a mesh with "solution" vectors, where the solution 
vector contains the data from the <code>MeshData</code>.  Particularly
useful for <i>inspecting</i> the data contained in <code> MeshData</code>.

<br><br>And even better: the user can choose the mesh for which
to export the data.  E.g. not only use the <code> mesh</code>
itself, but alternatively use the <code> BoundaryMesh </code>
(any mesh that uses the <i>same nodes</i>, i.e. the <code> Node*</code> 
have to be the same.  Only exception that will not work:
A mesh created using <code> Mesh::create_submesh() </code>
actually will <i>not</i> work with <code> MeshData::translate() </code>).

<br><br>All in all not bad, hm?
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
have a vector for the actual values and a vector
for the names of the data available.
</div>

<div class ="fragment">
<pre>
              std::vector&lt;Number&gt; translated_data;
              std::vector&lt;std::string&gt; data_names;
              
</pre>
</div>
<div class = "comment">
Use the <code> mesh</code> itself.  Alternatively, use the 
<code> BoundaryMesh</code> of <code> mesh</code>.
</div>

<div class ="fragment">
<pre>
              mesh_data.translate (mesh,
                                   translated_data,
                                   data_names);
        
        
</pre>
</div>
<div class = "comment">
And write the data to a GMV file
</div>

<div class ="fragment">
<pre>
              const std::string gmv_file = "data_and_mesh_out.gmv";
              std::cout &lt;&lt; std::endl 
                        &lt;&lt; "Writing the data from the MeshData to the GMV file " 
                        &lt;&lt; gmv_file &lt;&lt; std::endl
                        &lt;&lt; "------------------------------------------------------------------------" 
                        &lt;&lt; std::endl;
              
              GMVIO(mesh).write_nodal_data (gmv_file,
                                            translated_data,
                                            data_names);
              
</pre>
</div>
<div class = "comment">
Let the vectors with translated data
go out of scope.
</div>

<div class ="fragment">
<pre>
            }
            
</pre>
</div>
<div class = "comment">
Let the second mesh go out of scope.
</div>

<div class ="fragment">
<pre>
          }
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return 0;
        }
        
        
        
        
        
        
        
        
</pre>
</div>
<div class = "comment">
This function creates the data to populate the <code> MeshData</code> object
</div>

<div class ="fragment">
<pre>
        void create_artificial_data (const Mesh& mesh,
                                     std::map&lt;const Node*, std::vector&lt;Number&gt; &gt;& art_data)
        {
</pre>
</div>
<div class = "comment">
get the bounding box to have some sensible data
</div>

<div class ="fragment">
<pre>
          MeshTools::BoundingBox b_box = MeshTools::bounding_box(mesh);
        
          const Real z_min = b_box.first (2);
          const Real z_max = b_box.second(2);
          libmesh_assert (fabs(z_max-z_min) &gt; TOLERANCE);
        
          const Real x_min = b_box.first (0);
          const Real x_max = b_box.second(0);
          libmesh_assert (fabs(x_max-x_min) &gt; TOLERANCE);
        
        
</pre>
</div>
<div class = "comment">
const_node_iterator node_it = mesh.nodes_begin();
const const_node_iterator node_end = mesh.nodes_end();


<br><br></div>

<div class ="fragment">
<pre>
          MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
          const MeshBase::const_node_iterator node_end = mesh.nodes_end();
        
          for (; node_it != node_end; ++node_it)
            {
</pre>
</div>
<div class = "comment">
All the vectors in <code> artificial</code>_data <i>have</i> to have the
same size.  Here we use only two entries per node,
but theoretically arbitrary size is possible.
</div>

<div class ="fragment">
<pre>
              std::vector&lt;Number&gt; node_data;
              node_data.resize(2);
        
</pre>
</div>
<div class = "comment">
Use a sin curve in z-direction and a linear
increase in x-direction
</div>

<div class ="fragment">
<pre>
              const Point& p = **node_it;
                
              const Real z_normalized = (p(2)-z_min)/(z_max-z_min);
              const Real x_normalized = (p(0)-x_min)/(x_max-x_min);
        
              node_data[0] = sin(2*libMesh::pi*z_normalized);
              node_data[1] = x_normalized;
              
</pre>
</div>
<div class = "comment">
Insert the current data together with a pointer to
the current node in the map.
</div>

<div class ="fragment">
<pre>
              art_data.insert (std::make_pair(*node_it,node_data));
            }
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_tools.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_data.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;unv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;node.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> create_artificial_data (<B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh,
                               <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">const</FONT></B> Node*, std::vector&lt;Number&gt; &gt;&amp; art_data);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::n_processors() &gt; 1)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Skipping example 12. &quot;</FONT></B> &lt;&lt; std::endl;
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;MeshData objects currently only work in serial.&quot;</FONT></B> &lt;&lt; std::endl;
          }
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
      }
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 4)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -d &lt;dim&gt; in_mesh.unv&quot;</FONT></B>
                    &lt;&lt; std::endl;
        
        libmesh_error();
      }
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = std::atoi(argv[2]);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::string mesh_file = argv[3];
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (mesh_file.rfind(<B><FONT COLOR="#BC8F8F">&quot;.unv&quot;</FONT></B>) &gt;= mesh_file.size())
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR:  This example works only properly with a Universal mesh file!&quot;</FONT></B>
                  &lt;&lt; std::endl;
        
        libmesh_error();
      }
  
    {
      libmesh_assert (dim == 3);
      Mesh mesh;
      MeshData mesh_data(mesh);
    
      mesh_data.activate();
  
      mesh.read (mesh_file, &amp;mesh_data);
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Finished reading the mesh.  MeshData is active but empty:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;---------------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
      
      mesh.print_info();
      mesh_data.print_info();
    
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">const</FONT></B> Node*, std::vector&lt;Number&gt; &gt; artificial_data;
        
        create_artificial_data (mesh, artificial_data);
        
        mesh_data.insert_node_data(artificial_data);
  
      }
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;After inserting artificial data into the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;--------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
      
      mesh_data.print_info();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string first_out_data=<B><FONT COLOR="#BC8F8F">&quot;data_first_out_with_default_header.unv&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; first_out_data &lt;&lt; std::endl;
      mesh_data.write(first_out_data);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Attach our own MeshDataUnvHeader to the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;-------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; (note the warning: the number of values per node in&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  my_header is not correct)&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl;
      MeshDataUnvHeader my_header;
  
      my_header.dataset_label = 3;
  
      my_header.id_lines_1_to_5[0] = <B><FONT COLOR="#BC8F8F">&quot;Artificial data&quot;</FONT></B>;
      my_header.id_lines_1_to_5[1] = <B><FONT COLOR="#BC8F8F">&quot;sin curve in z-direction&quot;</FONT></B>;
      my_header.id_lines_1_to_5[2] = <B><FONT COLOR="#BC8F8F">&quot;line in x-direction&quot;</FONT></B>;
      
      my_header.record_12[0] = libMesh::pi;
      
      mesh_data.set_unv_header(&amp;my_header);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string second_out_data=<B><FONT COLOR="#BC8F8F">&quot;data_second_with_header_out.unv&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; second_out_data &lt;&lt; std::endl;
      mesh_data.write(second_out_data);
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Before clearing the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;-----------------------------&quot;</FONT></B> &lt;&lt; std::endl;
      mesh_data.print_info();
  
      mesh_data.clear();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;After clearing the MeshData:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;----------------------------&quot;</FONT></B> &lt;&lt; std::endl;
      mesh_data.print_info();
  
      mesh_data.read(first_out_data);
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;After re-reading the first file:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;--------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
      mesh_data.print_info();
  
    }
  
  
      
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;----------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;---------- next example with MeshData --------&quot;</FONT></B> &lt;&lt; std::endl
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;----------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl;
  
    {
      Mesh mesh(dim);
      MeshData mesh_data(mesh);
      
      mesh.read(mesh_file, &amp;mesh_data);
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;De-activated MeshData:&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;----------------------&quot;</FONT></B> &lt;&lt; std::endl;
      mesh.print_info();
      mesh_data.print_info();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::string out_mesh = <B><FONT COLOR="#BC8F8F">&quot;mesh_with_libmesh_ids.unv&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing _Mesh_ to: &quot;</FONT></B> &lt;&lt; out_mesh &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Try 'diff &quot;</FONT></B> &lt;&lt; out_mesh &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; mesh_file &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;'&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;to see the differences in node numbers.&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;---------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; std::endl;
      mesh.write(out_mesh, &amp;mesh_data);
  
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">const</FONT></B> Node*, std::vector&lt;Number&gt; &gt; artificial_data;
        create_artificial_data (mesh, artificial_data);
        mesh_data.insert_node_data(artificial_data);
      }
  
      mesh_data.enable_compatibility_mode();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string mesh_data_file = <B><FONT COLOR="#BC8F8F">&quot;data_third_with_libmesh_ids_out.unv&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing MeshData to: &quot;</FONT></B> &lt;&lt; mesh_data_file &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;----------------------------------------------------------&quot;</FONT></B> 
                &lt;&lt; std::endl &lt;&lt; std::endl;
      mesh_data.write (mesh_data_file);
  
  #ifdef HAVE_ZLIB_H
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string packed_mesh_data_file =  <B><FONT COLOR="#BC8F8F">&quot;packed_&quot;</FONT></B> + mesh_data_file + <B><FONT COLOR="#BC8F8F">&quot;.gz&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing gzip'ed MeshData to: &quot;</FONT></B> &lt;&lt; packed_mesh_data_file &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;---------------------------------------------------------------------------&quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; To verify the integrity of the packed version, type:&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;   gunzip &quot;</FONT></B> &lt;&lt; packed_mesh_data_file &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;; &quot;</FONT></B> &lt;&lt; std::endl
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;   diff packed_&quot;</FONT></B> &lt;&lt; mesh_data_file &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> 
                &lt;&lt; mesh_data_file &lt;&lt; std::endl &lt;&lt; std::endl;
      
      mesh_data.write (packed_mesh_data_file);
      
  #endif
  
  
        
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; translated_data;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; data_names;
        
        mesh_data.translate (mesh,
                             translated_data,
                             data_names);
  
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::string gmv_file = <B><FONT COLOR="#BC8F8F">&quot;data_and_mesh_out.gmv&quot;</FONT></B>;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl 
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing the data from the MeshData to the GMV file &quot;</FONT></B> 
                  &lt;&lt; gmv_file &lt;&lt; std::endl
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;------------------------------------------------------------------------&quot;</FONT></B> 
                  &lt;&lt; std::endl;
        
        GMVIO(mesh).write_nodal_data (gmv_file,
                                      translated_data,
                                      data_names);
        
      }
      
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> create_artificial_data (<B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh,
                               <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">const</FONT></B> Node*, std::vector&lt;Number&gt; &gt;&amp; art_data)
  {
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::BoundingBox b_box = MeshTools::bounding_box(mesh);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real z_min = b_box.first (2);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z_max = b_box.second(2);
    libmesh_assert (fabs(z_max-z_min) &gt; TOLERANCE);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real x_min = b_box.first (0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real x_max = b_box.second(0);
    libmesh_assert (fabs(x_max-x_min) &gt; TOLERANCE);
  
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator       node_it  = mesh.nodes_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_node_iterator node_end = mesh.nodes_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (; node_it != node_end; ++node_it)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; node_data;
        node_data.resize(2);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p = **node_it;
          
        <B><FONT COLOR="#228B22">const</FONT></B> Real z_normalized = (p(2)-z_min)/(z_max-z_min);
        <B><FONT COLOR="#228B22">const</FONT></B> Real x_normalized = (p(0)-x_min)/(x_max-x_min);
  
        node_data[0] = sin(2*libMesh::pi*z_normalized);
        node_data[1] = x_normalized;
        
        art_data.insert (std::make_pair(*node_it,node_data));
      }
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) ex12.C...
Linking ex12-opt...
***************************************************************
* Running Example  mpirun -np 2 ./ex12-opt -d 3 /h2/roystgnr/libmesh/svn/examples/ex8/pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
ERROR: Skipping example 12. 
MeshData objects currently only work in serial.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex12-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Tue Feb 22 12:19:38 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           5.019e-04      1.13845   4.714e-04
Objects:              0.000e+00      0.00000   0.000e+00
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Message Lengths:  0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Reductions:       0.000e+00      0.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.4537e-04  94.5%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.19209e-06
Average time for zero size MPI_Send(): 6.4373e-06
#PETSc Option Table entries:
-d 3
-ksp_right_pc
-log_summary
-pc_type bjacobi
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------
 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex12-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
</pre>
</div>
<?php make_footer() ?>
</body>
</html>
<?php if (0) { ?>
\#Local Variables:
\#mode: html
\#End:
<?php } ?>
