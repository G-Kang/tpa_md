#! /usr/bin/python

# This script converts a hoomd xml file to a LAMMPS input data file
import xml.dom.minidom
import sys

def writeout(path,infile,outfile):

    dom = xml.dom.minidom.parse(path+'/'+infile);

    # start by parsing the file
    hoomd_xml = dom.getElementsByTagName('hoomd_xml')[0];
    configuration = hoomd_xml.getElementsByTagName('configuration')[0];

    # read the box size
    box     = configuration.getElementsByTagName('box')[0];
    Lx = box.getAttribute('lx');
    Ly = box.getAttribute('ly');
    Lz = box.getAttribute('lz');

    # parse the particle coordinates
    position = configuration.getElementsByTagName('position')[0];
    position_text = position.childNodes[0].data
    xyz = position_text.split()
    print( "Found", len(xyz)//3, " particles" );

    # parse the particle types
    type_nodes = configuration.getElementsByTagName('type');
    if len(type_nodes) == 1:
        type_text = type_nodes[0].childNodes[0].data;
        type_names = type_text.split();
        if len(type_names) != len(xyz)//3:
            print( "Error! Number of types differes from the number of particles" )
            sys.exit(1);
    else:
        print( "Error! The type node must be in the xml file" )
        sys.exit(1);

    # parse the particle charges
    charge_nodes = configuration.getElementsByTagName('charge');
    if len(charge_nodes) == 1:
        charge_text = charge_nodes[0].childNodes[0].data;
        charge_names = charge_text.split();
        if len(charge_names) != len(xyz)//3:
            print( "Error! Number of charges differes from the number of particles" )
            sys.exit(1);
    else:
        print( "Error! The charge node must be in the xml file" )
        sys.exit(1);

    # parse the bonds
    bond_nodes = configuration.getElementsByTagName('bond')
    bond_a = [];
    bond_b = [];
    bond_type_id = [];
    bond_type_id_mapping = {};
    if len(bond_nodes) == 1:
        bond = bond_nodes[0];
        bond_text = bond.childNodes[0].data.encode();
        bond_raw = bond_text.split();
        
        # loop through the bonds and split the a,b and type from the raw stream
        # map types names to numbers along the way
        for i in range(0,len(bond_raw),3):
            bond_a.append(bond_raw[i+1]);
            bond_b.append(bond_raw[i+2]);
            
                
        print( "Found", len(bond_a), "bonds" )
        
    # now we have everything and can start writing the LAMMPS output file
    f = open(path+'/'+outfile, 'w');
    f.write("@<TRIPOS>MOLECULE\n%s\n" % infile.split('.')[0]);
    f.write("%d %d %d %d %d\nSMALL\nGASTEIGER\n" % (len(xyz)//3, len(bond_a),0,0,0));
    f.write("\n");
    f.write("@<TRIPOS>ATOM\n");
    for i in range(0,len(xyz)//3):
        f.write("%d %s %f %f %f %s %d %s %s\n" % (i+1, type_names[i],float(xyz[i*3]), float(xyz[i*3+1]), float(xyz[i*3+2]),type_names[i],1,"LIG1",charge_names[i]));
    if len(bond_a) > 0:
        f.write("@<TRIPOS>BOND\n");
        for i in range(0,len(bond_a)):
            f.write("%d %d %d %d\n" % (i+1, int(bond_a[i])+1, int(bond_b[i])+1,1));
    f.close()

