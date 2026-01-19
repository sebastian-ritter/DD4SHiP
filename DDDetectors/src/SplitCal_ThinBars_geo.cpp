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
  double       tol     = 80 * dd4hep::mm;
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
  const double thinbar_x_spacing   =  x_thinbar.attr<double>("x_extra_spacing");
  
  // Per-layer x_offsets: parse comma-separated list, or fall back to single x_offset
  std::vector<double> x_offsets;
  if (x_thinbar.hasAttr(_Unicode(x_offsets))) {
      std::string offsets_str = x_thinbar.attr<std::string>(_Unicode(x_offsets));
      x_offsets = parseOffsetList(offsets_str);
      printout(INFO, "SplitCal ThinBars", "%s: Parsed %zu per-layer x_offsets", nam.c_str(), x_offsets.size());
  } else {
      // Backward compatibility: use single x_offset for all layers
      double single_offset = x_thinbar.attr<double>("x_offset");
      x_offsets.push_back(single_offset);
      printout(INFO, "SplitCal ThinBars", "%s: Using single x_offset for all layers: %7.3f", nam.c_str(), single_offset);
  }
  const double y_offset   =  x_thinbar.attr<double>("y_offset");
  
  // Per-layer extrazgaps for passive layers: parse from passive_layer element
  std::vector<double> extrazgaps;
  if (x_passive_layer.hasAttr(_Unicode(extrazgaps))) {
      std::string gaps_str = x_passive_layer.attr<std::string>(_Unicode(extrazgaps));
      extrazgaps = parseOffsetList(gaps_str);
      printout(INFO, "SplitCal ThinBars", "%s: Parsed %zu per-layer passive extrazgaps", nam.c_str(), extrazgaps.size());
  } else {
      // Backward compatibility: use single extrazgap for all layers
      double single_gap = x_widebar.attr<double>("extrazgap");
      extrazgaps.push_back(single_gap);
      printout(INFO, "SplitCal ThinBars", "%s: Using single passive extrazgap for all layers: %7.3f", nam.c_str(), single_gap);
  }
  
  // Per-layer extrazgaps for active thin bar layers
  std::vector<double> thinbar_extrazgaps;
  if (x_thinbar.hasAttr(_Unicode(extrazgaps))) {
      std::string gaps_str = x_thinbar.attr<std::string>(_Unicode(extrazgaps));
      thinbar_extrazgaps = parseOffsetList(gaps_str);
      printout(INFO, "SplitCal ThinBars", "%s: Parsed %zu per-layer thinbar extrazgaps", nam.c_str(), thinbar_extrazgaps.size());
  } else {
      // Backward compatibility: use single extrazgap for all layers
      double single_gap = x_thinbar.attr<double>("extrazgap");
      thinbar_extrazgaps.push_back(single_gap);
      printout(INFO, "SplitCal ThinBars", "%s: Using single thinbar extrazgap for all layers: %7.3f", nam.c_str(), single_gap);
  }
  const std::string calo_layer_codes = x_det.attr<std::string>("layer_codes");
  const int num_z   =  static_cast<unsigned>(calo_layer_codes.size()); 
  const int thinbar_num_x   =  x_thinbar.attr<unsigned>("num_x");

  //HPL fibre feature extraction 
  xml_dim_t    x_hplbox   = x_det.child(_Unicode(hplbox));
  
  
  //Bar definition
  Box   thinbar((x_thinbar.x()-tol)/2., (x_thinbar.y()-tol)/2.,(x_thinbar.z()-tol)/2.);
  Box   passive_layer_box((x_passive_layer.x()-tol)/2., (x_passive_layer.y()-tol)/2.,(x_passive_layer.z()-tol)/2.);
  Box   split_box((x_split.x()-tol)/2., (x_split.y()-tol)/2.,(x_split.z()-tol)/2.);
  Volume thinbar_vol("thinbar", thinbar, description.material(x_thinbar.materialStr()));
  Volume passive_layer_vol("passive_layer", passive_layer_box, description.material(x_passive_layer.materialStr()));
  Volume split_vol("split", split_box, description.material(x_split.materialStr()));
  thinbar_vol.setAttributes(description, x_thinbar.regionStr(), x_thinbar.limitsStr(), x_thinbar.visStr());
  passive_layer_vol.setAttributes(description, x_passive_layer.regionStr(), x_passive_layer.limitsStr(), x_passive_layer.visStr());
  split_vol.setAttributes(description, x_split.regionStr(), x_split.limitsStr(), x_split.visStr());


  sens.setType("calorimeter");

  // Envelope: make envelope box 'tol' bigger on each side
  
  Assembly detbox_vol(nam + "_assembly");
  
  //Box    detbox((x_detbox.x()+tol)/2., (x_detbox.y()+tol)/2., (x_detbox.z()+tol)/2.);
  //Volume detbox_vol(nam, detbox, description.air());
  //detbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  
//  box_vol.setVisAttributes(description.visAttributes(""));
 
  double thinlayerwidth = x_thinbar.x() * static_cast<double>(thinbar_num_x);

  Box    det_thin_layerbox((thinlayerwidth+tol)/2., (x_thinbar.y()+tol)/2., (x_thinbar.z()+tol)/2.);
  Volume det_thin_layerbox_vol("det_thin_layerbox", det_thin_layerbox, description.air());
  det_thin_layerbox_vol.setAttributes(description, x_detbox.regionStr(), x_detbox.limitsStr(), x_detbox.visStr());
  det_thin_layerbox_vol.setVisAttributes(description.visAttributes(x_detbox.visStr()));
  
//  Rotation3D rot(RotationZYX(0e0, 0e0, M_PI/2e0));
  Rotation3D rot(RotationZYX(0e0, 0e0, 0e0));
  
  //if( x_thinbar.hasChild(_U(sensitive)) )  {
  //  sens.setType("calorimeter");
    thinbar_vol.setSensitiveDetector(sens);
  //}
  
  //if(x_hplcore.hasChild(_U(sensitive)) )  {
  //  sens.setType("calorimeter");
  //}
  //Loop for x-wise placement -> build the sensitive bar layer 
  //

  //Volume encoding
  //long DetectorCode = 9 * 1e15; 
  //long ECALCode = 1 * 1e12;

//  int DetectorCode = 9 * 1e8; 
//  int ECALCode = 1 * 1e7;

  //Thin bar layers
  int volumecode = 0;
  double xpos = -(thinlayerwidth+tol)/2.;
  for( int ix=0; ix < thinbar_num_x; ++ix )  {
    xpos += x_thinbar.x()/2.; 
    PlacedVolume pv = det_thin_layerbox_vol.placeVolume(thinbar_vol, Transform3D(rot,Position(xpos, 0e0, 0e0)));
    pv.addPhysVolID("splitcal_thin_bar", volumecode);
    xpos += x_thinbar.x()/2. + thinbar_x_spacing; 
    volumecode++;
  }



//HPL Layers

  //Definition of layer volumes


 
  
  double z_layer = -x_detbox.z()/2.;
  Rotation3D rot_layers;
  int thin_layer_count = 0;  // Counter for thin bar layers (to index x_offsets)
  int passive_layer_count = 0;  // Counter for passive layers (to index extrazgaps)


  for( int iz=0; iz < num_z; ++iz )  {
    // leave 'tol' space between the layers

    
    std::cout << static_cast<int>(calo_layer_codes[iz]) - '0' << std::endl;
    switch(static_cast<int>(calo_layer_codes[iz]) - '0'){
    
	    case 1:{
		//Leave space for wide bars
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 2:{
		//Leave space for wide bars
    		z_layer += x_widebar.z();
   		break;
	    }
	    case 3:{
		//Place thin layer vertically
    		z_layer += x_thinbar.z()/2.;
	    	// rot_layers = RotationZYX(M_PI/2e0,0e0,0e0);
		// Get per-layer offset (cycle if fewer offsets than layers)
		double layer_x_offset = x_offsets[thin_layer_count % x_offsets.size()];
		// DEBUG: Print absolute x positions of each bar in this layer
		printout(INFO, "SplitCal ThinBars", "Layer %d (code 3): z=%7.2f mm, x_offset=%7.2f mm", thin_layer_count, z_layer, layer_x_offset);
		double bar_xpos = -(thinlayerwidth)/2.;
		for (int ib = 0; ib < thinbar_num_x; ++ib) {
		    bar_xpos += x_thinbar.x()/2.;
		    double abs_x = layer_x_offset + bar_xpos;
		    printout(INFO, "SplitCal ThinBars", "  Bar %d: local_x=%7.2f mm, absolute_x=%7.2f mm", ib, bar_xpos, abs_x);
		    bar_xpos += x_thinbar.x()/2.;
		}
    	    	PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(layer_x_offset,y_offset , z_layer)));
    	    	pv_det.addPhysVolID("splitcal_thin_layer", iz);
    		z_layer += x_thinbar.z()/2.;
		// Add extrazgap after active layer
		double layer_extrazgap_thin = thinbar_extrazgaps[thin_layer_count % thinbar_extrazgaps.size()];
		z_layer += layer_extrazgap_thin;
		thin_layer_count++;
   		break;
            }
	    case 4:{
		//Place thin layer horizontally	
    		z_layer += x_thinbar.z()/2.;
		rot_layers = RotationZYX(0e0, 0e0, 0e0);
		// Get per-layer offset (cycle if fewer offsets than layers)
		double layer_x_offset = x_offsets[thin_layer_count % x_offsets.size()];
    		PlacedVolume pv_det = detbox_vol.placeVolume(det_thin_layerbox_vol, Transform3D(rot_layers,Position(y_offset,layer_x_offset, z_layer)));
        	pv_det.addPhysVolID("splitcal_thin_layer", iz);
    		z_layer += x_thinbar.z()/2.;
		// Add extrazgap after active layer
		double layer_extrazgap_thin4 = thinbar_extrazgaps[thin_layer_count % thinbar_extrazgaps.size()];
		z_layer += layer_extrazgap_thin4;
		thin_layer_count++;
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
	    case 7:{
		//Leave space for passive layers
    		z_layer += x_passive_layer.z();
    		// Get per-layer extrazgap (cycle if fewer gaps than passive layers)
    		double layer_extrazgap = extrazgaps[passive_layer_count % extrazgaps.size()];
    		z_layer += layer_extrazgap;
    		passive_layer_count++;
		break;
		   }
	    case 8:{
		//Leave space for split 
    		z_layer += x_split.z();
		break;
		   }
    }
    std::cout << "Zlayer Det " << z_layer << std::endl;

  }
  
  DetElement   sdet  (nam, x_det.id());
  Volume       mother(description.pickMotherVolume(sdet));
  Rotation3D   rot3D (RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)));
  Transform3D  trafo (rot3D, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
  PlacedVolume pv = mother.placeVolume(detbox_vol, trafo);
  pv.addPhysVolID("system", x_det.id());
  sdet.setPlacement(pv);  // associate the placed volume to the detector element
  printout(INFO, "SplitCal Thinbars", "%s: Detector construction finished.", nam.c_str());
  return sdet;
}

DECLARE_DETELEMENT(DD4hep_SplitCalThinBars,create_detector)
