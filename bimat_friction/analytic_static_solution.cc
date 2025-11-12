#include <iostream>
#include <cmath>

#include "solid_mechanics_model.hh"
#include "tasn_contact.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  debug::setDebugLevel(dblWarning);
  initialize(argv[1], argc, argv);

  // Communicator initialization
  //StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  //Int psize = comm.getNbProc();
  Int prank = 0; //comm.whoAmI();

  if (prank == 0)
    std::cout << "simulation_code = akantu" << std::endl;

  const ParserSection & data = getUserParser();
  std::string simulation_name = data.getParameter("simulation_name");

  std::stringstream output_folder;
  output_folder << data.getParameterValue<std::string>("output_folder");
  if (prank == 0)
    std::cout << "output_folder = " << output_folder.str() << std::endl;

  std::stringstream paraview_folder;
  paraview_folder << data.getParameterValue<std::string>("paraview_folder");
  if (prank == 0)
    std::cout << "paraview_folder = " << paraview_folder.str() << std::endl;

  std::stringstream restart_dump_folder;
  restart_dump_folder << data.getParameterValue<std::string>("restart_dump_folder");
  if (prank == 0)
    std::cout << "restart_dump_folder = " << restart_dump_folder.str() << std::endl;

  // copy input file to output folder
  std::ifstream src_ifile(getStaticParser().getLastParsedFile(), std::ios::binary);
  std::ofstream dst_ifile(output_folder.str() + simulation_name + ".in", std::ios::binary);
  dst_ifile << src_ifile.rdbuf();
  src_ifile.close();
  dst_ifile.close();

  UInt spatial_dimension = data.getParameter("spatial_dimension");

  // MESH & GEOMETRY
  Mesh mesh(spatial_dimension);
  mesh.read(data.getParameter("mesh"));
  mesh.computeBoundingBox();

  // ----------------------------------------------------------------------------------------------------

  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);
  model.initArrays();

  // get all materials sections from the parser
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sects = getStaticParser().getSubSections(_st_material);

  // loop over all element groups and see if a material section with corresponding name exists
  for(Mesh::const_element_group_iterator it(mesh.element_group_begin()); 
      it != mesh.element_group_end(); ++it) {
    const ElementGroup eg = *(it->second);

    Parser::const_section_iterator sit = sub_sects.first;
    for (; sit != sub_sects.second; ++sit) {
      if (sit->getParameterValue<std::string>("name") == eg.getName()) {
	std::cout << "Found ElementGroup '" << eg.getName() << "' with material" << std::endl;
	break;
      }
    }

    // found a match for element group and material
    // apply equilibrium state to the nodes of this element group
    if (sit != sub_sects.second) {
      const ParserSection & mat = *sit;
      Real E  = mat.getParameter("E");
      Real nu = mat.getParameter("nu");

      if (spatial_dimension == 2) {
	// adapt E and nu for plane stress computation
	if (mat.hasParameter("Plane_Stress") 
	    && mat.getParameter("Plane_Stress")) {
	  std::cout << "Plane stress!! "  << std::endl;
	  E *= (1.+2.*nu)/((1.+nu)*(1.+nu));
	  nu = nu/(1.0+nu);
	}
	else
	  std::cout << "Plane strain!! "  << std::endl;

	Vector<Real> trac = data.getParameter("top_traction");

	Real shear_modulus = E / (2. * (1. + nu));
	Real shear_strain_applied = trac(0) / (2*shear_modulus);
	std::cout << "shear strain = " << shear_strain_applied << std::endl;

	Real normal_strain_applied = trac(1) * (1+nu)*(1-nu) / E;
	Real normal_strain_applied_2 = trac(1) * (-nu) * (1+nu) / E;

	Array<Real> & disp = model.getDisplacement();
	const Array<Real> & pos = mesh.getNodes();
	const Array<UInt> & group_nodes = eg.getNodes();

	for (UInt i = 0; i < group_nodes.getSize(); ++i) {
	  Real xp = pos(group_nodes(i),0);
	  Real yp = pos(group_nodes(i),1);
	  disp(group_nodes(i), 0) = shear_strain_applied * yp + normal_strain_applied_2 * xp;
	  disp(group_nodes(i), 1) = shear_strain_applied * xp + normal_strain_applied   * yp;
	}
      }
      else if (spatial_dimension == 3) {
	
	Vector<Real> trac = data.getParameter("top_traction");
	
	Real shear_modulus = E / (2. * (1. + nu));
	Real shear_strain_applied = trac(0) / (2*shear_modulus);
	std::cout << "shear strain = " << shear_strain_applied << std::endl;
	
	Real normal_strain_applied = trac(1) / E;
	Real normal_strain_applied_2 = trac(1) * (-nu) / E;
	
	Array<Real> & disp = model.getDisplacement();
	const Array<Real> & pos = mesh.getNodes();
	const Array<UInt> & group_nodes = eg.getNodes();
	
	for (UInt i = 0; i < group_nodes.getSize(); ++i) {
	  Real xp = pos(group_nodes(i),0);
	  Real yp = pos(group_nodes(i),1);
	  Real zp = pos(group_nodes(i),2);
	  disp(group_nodes(i), 0) = shear_strain_applied * yp + normal_strain_applied_2 * xp;
	  disp(group_nodes(i), 1) = shear_strain_applied * xp + normal_strain_applied   * yp;
	  disp(group_nodes(i), 2) =                             normal_strain_applied_2 * zp;
	}
      }
    }
  }

  // dump restart file
  dumpArray(model.getDisplacement(),
	    restart_dump_folder.str() 
	    + simulation_name + ".restart");

  std::cout << "analytic_static_solution went until the end" << std::endl;
  finalize();
  return EXIT_SUCCESS;
}
