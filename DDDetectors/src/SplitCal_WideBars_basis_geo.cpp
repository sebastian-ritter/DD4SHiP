//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Climescu
// Date       : 09.11.2025
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Printout.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace dd4hep;

// Helper function to parse a comma-separated list of offset values with units
// Example input: "0*cm, 1*cm, -0.5*cm, 2*mm"
// Returns a vector of doubles in DD4hep internal units
static std::vector<double> parseOffsetList(const std::string& input) {
    std::vector<double> result;
    std::stringstream ss(input);
    std::string token;
    
    while (std::getline(ss, token, ',')) {
        // Trim whitespace
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start == std::string::npos) continue;
        token = token.substr(start, end - start + 1);
        
        // Parse value and unit (e.g., "1.5*cm" or "-0.5*mm")
        size_t mult_pos = token.find('*');
        if (mult_pos != std::string::npos) {
            double value = std::stod(token.substr(0, mult_pos));
            std::string unit = token.substr(mult_pos + 1);
            
            // Convert to DD4hep internal units
            if (unit == "cm") value *= dd4hep::cm;
            else if (unit == "mm") value *= dd4hep::mm;
            else if (unit == "m") value *= dd4hep::m;
            // If no known unit, assume it's already in internal units
            
            result.push_back(value);
        } else {
            // No unit specified, assume internal units
            result.push_back(std::stod(token));
        }
    }
    return result;
}

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  //Calo scintillator bars' feature extraction
  double       tol     = 0 * dd4hep::mm;
  xml_det_t    x_det   = e;
  xml_dim_t    x_detbox   = x_det.child(_U(box));
  xml_dim_t    x_rot   = x_det.child(_U(rotation));
  xml_dim_t    x_pos   = x_det.child(_U(position));
  xml_det_t    x_widebar = x_det.child(_Unicode(widebar));
  xml_det_t    x_thinbar = x_det.child(_Unicode(thinbar));
  xml_det_t    x_passive_layer = x_det.child(_Unicode(passive_layer));
  xml_det_t    x_split = x_det.child(_Unicode(split));
  std::string  nam     = x_det.nameStr();
  //vertical bars by default
//  const double splitlayer   =  x_det.attr<int>("splitlayer");
  
  // Per-layer x_offsets: parse comma-separated list, or fall back to single x_offset
  std::vector<double> x_offsets;
  if (x_widebar.hasAttr(_Unicode(x_offsets))) {
      std::string offsets_str = x_widebar.attr<std::string>(_Unicode(x_offsets));
      x_offsets = parseOffsetList(offsets_str);
      printout(INFO, "SplitCal", "%s: Parsed %zu per-layer x_offsets", nam.c_str(), x_offsets.size());
  } else {
      // Backward compatibility: use single x_offset for all layers
      double single_offset = x_widebar.attr<double>("x_offset");
      x_offsets.push_back(single_offset);
      printout(INFO, "SplitCal", "%s: Using single x_offset for all layers: %7.3f", nam.c_str(), single_offset);
  }
  const double y_offset = x_widebar.attr<double>("y_offset");
  const double widebar_x_spacing   =  x_widebar.attr<double>("x_extra_spacing");
  
  // Per-layer extrazgaps for passive layers: parse from passive_layer element
  std::vector<double> extrazgaps;
  if (x_passive_layer.hasAttr(_Unicode(extrazgaps))) {
      std::string gaps_str = x_passive_layer.attr<std::string>(_Unicode(extrazgaps));
      extrazgaps = parseOffsetList(gaps_str);
      printout(INFO, "SplitCal", "%s: Parsed %zu per-layer passive extrazgaps", nam.c_str(), extrazgaps.size());
  } else {
      // Backward compatibility: use single extrazgap for all layers
      double single_gap = x_widebar.attr<double>("extrazgap");
      extrazgaps.push_back(single_gap);
      printout(INFO, "SplitCal", "%s: Using single passive extrazgap for all layers: %7.3f", nam.c_str(), single_gap);
  }
  
  // Per-layer extrazgaps for active wide bar layers
  std::vector<double> widebar_extrazgaps;
  if (x_widebar.hasAttr(_Unicode(extrazgaps))) {
      std::string gaps_str = x_widebar.attr<std::string>(_Unicode(extrazgaps));
      widebar_extrazgaps = parseOffsetList(gaps_str);
      printout(INFO, "SplitCal", "%s: Parsed %zu per-layer widebar extrazgaps", nam.c_str(), widebar_extrazgaps.size());
  } else {
      // Backward compatibility: use single extrazgap for all layers
      std::string gap_str = x_widebar.attr<std::string>("extrazgap");
      printout(INFO, "SplitCal", "%s: DEBUG - Raw extrazgap string: '%s' (length=%zu)", nam.c_str(), gap_str.c_str(), gap_str.length());
      std::vector<double> parsed_gap = parseOffsetList(gap_str);
      printout(INFO, "SplitCal", "%s: DEBUG - Parsed gap vector size: %zu", nam.c_str(), parsed_gap.size());
      if (!parsed_gap.empty()) {
          widebar_extrazgaps.push_back(parsed_gap[0]);
          printout(INFO, "SplitCal", "%s: DEBUG - Parsed value (internal units): %7.3f, dd4hep::cm=%7.3f, dd4hep::mm=%7.3f", nam.c_str(), parsed_gap[0], dd4hep::cm, dd4hep::mm);
      }
      printout(INFO, "SplitCal", "%s: Using single widebar extrazgap for all layers: %7.3f", nam.c_str(), widebar_extrazgaps.empty() ? 0.0 : widebar_extrazgaps[0]);
    }
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int widebar_num_x   =  x_widebar.attr<unsigned>("num_x");

  //HPL fibre feature extraction 
  xml_dim_t    x_hplbox   = x_det.child(_Unicode(hplbox));
//  const int    num_z   = int(2e0*x_box.z() / (delta+2*tol));
  
  
  //Bar definition
  Box   widebar((x_widebar.x()-tol)/2., (x_widebar.y()-tol)/2.,(x_widebar.z()-tol)/2.);
  Box   passive_layer_box((x_passive_layer.x()-tol)/2., (x_passive_layer.y()-tol)/2.,(x_passive_layer.z()-tol)/2.);
  Box   split_box((x_split.x()-tol)/2., (x_split.y()-tol)/2.,(x_split.z()-tol)/2.);
  Volume widebar_vol("widebar", widebar, description.material(x_widebar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  Volume split_vol("split", split_box, description.material(x_split.materialStr()));
  widebar_vol.setAttributes(description, x_widebar.regionStr(), x_widebar.limitsStr(), x_widebar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());
  split_vol.setAttributes(description, x_split.regionStr(), x_split.limitsStr(), x_split.visStr());

  printout(INFO, "SandwichCalo", "%s: Bars: x: %7.3f y: %7.3f z: %7.3f mat: %s vis: %s solid: %s",
           nam.c_str(), x_widebar.x(), x_widebar.y(), x_widebar.z(), x_widebar.materialStr().c_str(),
           x_widebar.visStr().c_str(), widebar.type());
  sens.setType("calorimeter");

  // Envelope: make envelope box 'tol' bigger on each side
  //Box    detbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_detbox.z()+tol)/2.);
  //Volume detbox_vol(nam, detbox, description.air());
  //detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
 
  Assembly detbox_vol(nam + "_assembly");	

//  box_vol.setVisAttributes(description.visAttributes(""));
	
  double wideboxwidth = x_widebar.x()*static_cast<double>(widebar_num_x);
  Box    det_wide_layerbox((wideboxwidth+tol)/2., (x_widebar.y()+tol)/2., (x_widebar.z()+tol)/2.);
  Volume det_wide_layerbox_vol("det_wide_layerbox", det_wide_layerbox, description.air());
  det_wide_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_wide_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
  printout(INFO, "SandwichCalo", "%s: Layer:   nx: %7d x spacing: %7.3f", nam.c_str(), widebar_num_x, widebar_x_spacing);
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
 // if( x_widebar.hasChild(_U(sensitive)) )  {
  //  sens.setType("calorimeter");
    widebar_vol.setSensitiveDetector(sens);
  //}

  //Build Wide bar layers
  double xpos = -(wideboxwidth+tol)/2.;
  int volumecode = 0;
  for( int ix=0; ix < widebar_num_x; ++ix )  {
 
    xpos += x_widebar.x()/2.;
    PlacedVolume pv = det_wide_layerbox_vol.placeVolume(widebar_vol, Transform3D(rot,Position(xpos, 0e0, 0e0)));
    pv.addPhysVolID("splitcal_wide_bar", volumecode);
    xpos += x_widebar.x()/2. +widebar_x_spacing; 
    volumecode++;
  }
  xpos = -x_detbox.x()/2.;



//HPL Layers

  //Definition of layer volumes

  //Loop for z-wide placement -> build the calorimeter sandwich
  
  double z_layer = -x_detbox.z()/2.;
  Rotation3D rot_layers;
  int wide_layer_count = 0;  // Counter for wide bar layers (to index x_offsets)
  int passive_layer_count = 0;  // Counter for passive layers (to index extrazgaps)


  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    //std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Place wide layer vertically
    		z_layer += x_widebar.z()/2.;
	    	rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
		// Get per-layer offset (cycle if fewer offsets than layers)
		double layer_x_offset = x_offsets[wide_layer_count % x_offsets.size()];
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(layer_x_offset,y_offset,z_layer)));
    	    	pv_det.addPhysVolID("splitcal_wide_layer", iz);
    		z_layer += x_widebar.z()/2.;
		// Add extrazgap after active layer
		double layer_extrazgap_wide = widebar_extrazgaps[wide_layer_count % widebar_extrazgaps.size()];
		z_layer += layer_extrazgap_wide;
		wide_layer_count++;
   		break;
	    }
	    case 2:{
		//Place wide layer horizontally	
    		z_layer += x_widebar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
		// Get per-layer offset (cycle if fewer offsets than layers)
		double layer_x_offset = x_offsets[wide_layer_count % x_offsets.size()];
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_wide_layerbox_vol, Transform3D(rot_layers,Position(layer_x_offset,y_offset,z_layer)));
        	pv_det.addPhysVolID("splitcal_wide_layer", iz);
    		z_layer += x_widebar.z()/2.;
		// Add extrazgap after active layer
		double layer_extrazgap_wide2 = widebar_extrazgaps[wide_layer_count % widebar_extrazgaps.size()];
		z_layer += layer_extrazgap_wide2;
		wide_layer_count++;
	  	break; 
	    }
	    case 3:{
		//Leave space for thin bars
    		z_layer += x_thinbar.z();
   		break;
            }
	    case 4:{
		//Leave space for thin bars
    		z_layer += x_thinbar.z();
		break;
            }
	    case 5:{
		//Leave space for HPL
		z_layer += x_hplbox.z();
		break;
		   }
	    case 6:{
		//Leave space for HPL
		z_layer += x_hplbox.z();
		break;
		   }
	    case 7:{//Place passive layer
    		z_layer += x_passive_layer.z()/2.;
    		PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(0.,0., z_layer)));
    		pv_passive.addPhysVolID("splitcal_passivelayer", iz);
    		z_layer += x_passive_layer.z()/2.;
    		// Get per-layer extrazgap (cycle if fewer gaps than passive layers)
    		double layer_extrazgap = extrazgaps[passive_layer_count % extrazgaps.size()];
    		z_layer += layer_extrazgap;
    		passive_layer_count++;
		break;
		   }
	    case 8:{//Place split
    		z_layer += x_split.z()/2.;
    		PlacedVolume pv_split = detbox_vol.placeVolume(split_vol,Transform3D(rot,Position(0.,0. , z_layer)));
    		pv_split.addPhysVolID("splitcal_split_layer", iz);
    		z_layer += x_split.z()/2.;
		break;
		   }
    }
    std::cout << "Zlayer Det (WideBars) " << z_layer << std::endl;

//    if(static_cast<int>(calo_layer_codes[iz]) - '0' != 5 && static_cast<int>(calo_layer_codes[iz]) - '0' != 6){
//    	z_layer += x_passive_layer.z()/2.;
//    	PlacedVolume pv_passive = detbox_vol.placeVolume(passive_layer_vol,Transform3D(rot_layers,Position(x_passive_layer.x(), 0e0, z_layer)));
//    	z_layer += x_passive_layer.z()/2.;
//    	std::cout << "Zlayer passive " << z_layer << std::endl;
//    	pv_passive.addPhysVolID("passivelayer", iz*2+1);
//    	z_layer += extrazgap;
//    }
//    if(iz==splitlayer){
//    	z_layer += x_split.z()+x_widebar.z()+x_passive_layer.z();
//    	PlacedVolume pv_split = detbox_vol.placeVolume(split_vol,Transform3D(rot,Position(x_split.x(), 0e0, z_layer)));
//    	pv_split.addPhysVolID("split_layer", 9000);
//    }
  }
  std::cout << "Zlayer Det (WideBars) Final: " << z_layer << std::endl;

//  printout(INFO, "SandwichCalo", "%s: Created %d layers of %d bars each.", nam.c_str(), num_z, num_x);
  //PlacedVolume pv2 = detbox_vol.placeVolume(det_layerbox_vol, Transform3D(rot,Position(0e0, 0e0, 0e0)));
  //pv2.addPhysVolID("det_layerbox", 0e0);
  
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
// PlacedVolume pv2 = mother.placeVolume(passive_layer_vol, trafo);
  PlacedVolume pv = mother.placeVolume(detbox_vol, trafo);
  //pv2.addPhysVolID("system", x_det.id());
  pv.addPhysVolID("system", x_det.id());
  //sdet.setPlacement(pv2);  // associate the placed volume to the detector element
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "SplitCal", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SplitCalWideBars_and_Basis,create_detector)
