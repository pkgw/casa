
/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File CalAntennaSolutionsRow.h
 */
 
#ifndef CalAntennaSolutionsRow_CLASS
#define CalAntennaSolutionsRow_CLASS

#include <vector>
#include <string>
#include <set>

#ifndef WITHOUT_ACS
#include <asdmIDLC.h>
#endif






	 
#include <Angle.h>
	

	 
#include <ArrayTime.h>
	

	 
#include <Frequency.h>
	

	 
#include <Tag.h>
	

	 
#include <Interval.h>
	




	

	
#include "CAtmPhaseCorrection.h"
	

	
#include "CBasebandName.h"
	

	
#include "CReceiverBand.h"
	

	

	

	

	

	

	

	

	
#include "CPolarizationType.h"
	

	

	

	

	

	



#include <ConversionException.h>
#include <NoSuchRow.h>
#include <IllegalAccessException.h>

#include <RowTransformer.h>
//#include <TableStreamReader.h>

/*\file CalAntennaSolutions.h
    \brief Generated from model's revision "-1", branch ""
*/

namespace asdm {

//class asdm::CalAntennaSolutionsTable;


// class asdm::CalDataRow;
class CalDataRow;

// class asdm::CalReductionRow;
class CalReductionRow;
	

class CalAntennaSolutionsRow;
typedef void (CalAntennaSolutionsRow::*CalAntennaSolutionsAttributeFromBin) (EndianIStream& eis);
typedef void (CalAntennaSolutionsRow::*CalAntennaSolutionsAttributeFromText) (const string& s);

/**
 * The CalAntennaSolutionsRow class is a row of a CalAntennaSolutionsTable.
 * 
 * Generated from model's revision "-1", branch ""
 *
 */
class CalAntennaSolutionsRow {
friend class asdm::CalAntennaSolutionsTable;
friend class asdm::RowTransformer<CalAntennaSolutionsRow>;
//friend class asdm::TableStreamReader<CalAntennaSolutionsTable, CalAntennaSolutionsRow>;

public:

	virtual ~CalAntennaSolutionsRow();

	/**
	 * Return the table to which this row belongs.
	 */
	CalAntennaSolutionsTable &getTable() const;
	
	/**
	 * Has this row been added to its table ?
	 * @return true if and only if it has been added.
	 */
	bool isAdded() const;
		
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaName
	
	
	

	
 	/**
 	 * Get antennaName.
 	 * @return antennaName as string
 	 */
 	string getAntennaName() const;
	
 
 	
 	
 	/**
 	 * Set antennaName with the specified string.
 	 * @param antennaName The string value to which antennaName is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setAntennaName (string antennaName);
  		
	
	
	


	
	// ===> Attribute atmPhaseCorrection
	
	
	

	
 	/**
 	 * Get atmPhaseCorrection.
 	 * @return atmPhaseCorrection as AtmPhaseCorrectionMod::AtmPhaseCorrection
 	 */
 	AtmPhaseCorrectionMod::AtmPhaseCorrection getAtmPhaseCorrection() const;
	
 
 	
 	
 	/**
 	 * Set atmPhaseCorrection with the specified AtmPhaseCorrectionMod::AtmPhaseCorrection.
 	 * @param atmPhaseCorrection The AtmPhaseCorrectionMod::AtmPhaseCorrection value to which atmPhaseCorrection is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setAtmPhaseCorrection (AtmPhaseCorrectionMod::AtmPhaseCorrection atmPhaseCorrection);
  		
	
	
	


	
	// ===> Attribute basebandName
	
	
	

	
 	/**
 	 * Get basebandName.
 	 * @return basebandName as BasebandNameMod::BasebandName
 	 */
 	BasebandNameMod::BasebandName getBasebandName() const;
	
 
 	
 	
 	/**
 	 * Set basebandName with the specified BasebandNameMod::BasebandName.
 	 * @param basebandName The BasebandNameMod::BasebandName value to which basebandName is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setBasebandName (BasebandNameMod::BasebandName basebandName);
  		
	
	
	


	
	// ===> Attribute receiverBand
	
	
	

	
 	/**
 	 * Get receiverBand.
 	 * @return receiverBand as ReceiverBandMod::ReceiverBand
 	 */
 	ReceiverBandMod::ReceiverBand getReceiverBand() const;
	
 
 	
 	
 	/**
 	 * Set receiverBand with the specified ReceiverBandMod::ReceiverBand.
 	 * @param receiverBand The ReceiverBandMod::ReceiverBand value to which receiverBand is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setReceiverBand (ReceiverBandMod::ReceiverBand receiverBand);
  		
	
	
	


	
	// ===> Attribute startValidTime
	
	
	

	
 	/**
 	 * Get startValidTime.
 	 * @return startValidTime as ArrayTime
 	 */
 	ArrayTime getStartValidTime() const;
	
 
 	
 	
 	/**
 	 * Set startValidTime with the specified ArrayTime.
 	 * @param startValidTime The ArrayTime value to which startValidTime is to be set.
 	 
 		
 			
 	 */
 	void setStartValidTime (ArrayTime startValidTime);
  		
	
	
	


	
	// ===> Attribute endValidTime
	
	
	

	
 	/**
 	 * Get endValidTime.
 	 * @return endValidTime as ArrayTime
 	 */
 	ArrayTime getEndValidTime() const;
	
 
 	
 	
 	/**
 	 * Set endValidTime with the specified ArrayTime.
 	 * @param endValidTime The ArrayTime value to which endValidTime is to be set.
 	 
 		
 			
 	 */
 	void setEndValidTime (ArrayTime endValidTime);
  		
	
	
	


	
	// ===> Attribute numReceptor
	
	
	

	
 	/**
 	 * Get numReceptor.
 	 * @return numReceptor as int
 	 */
 	int getNumReceptor() const;
	
 
 	
 	
 	/**
 	 * Set numReceptor with the specified int.
 	 * @param numReceptor The int value to which numReceptor is to be set.
 	 
 		
 			
 	 */
 	void setNumReceptor (int numReceptor);
  		
	
	
	


	
	// ===> Attribute refAntennaName
	
	
	

	
 	/**
 	 * Get refAntennaName.
 	 * @return refAntennaName as string
 	 */
 	string getRefAntennaName() const;
	
 
 	
 	
 	/**
 	 * Set refAntennaName with the specified string.
 	 * @param refAntennaName The string value to which refAntennaName is to be set.
 	 
 		
 			
 	 */
 	void setRefAntennaName (string refAntennaName);
  		
	
	
	


	
	// ===> Attribute direction
	
	
	

	
 	/**
 	 * Get direction.
 	 * @return direction as vector<Angle >
 	 */
 	vector<Angle > getDirection() const;
	
 
 	
 	
 	/**
 	 * Set direction with the specified vector<Angle >.
 	 * @param direction The vector<Angle > value to which direction is to be set.
 	 
 		
 			
 	 */
 	void setDirection (vector<Angle > direction);
  		
	
	
	


	
	// ===> Attribute frequencyRange
	
	
	

	
 	/**
 	 * Get frequencyRange.
 	 * @return frequencyRange as vector<Frequency >
 	 */
 	vector<Frequency > getFrequencyRange() const;
	
 
 	
 	
 	/**
 	 * Set frequencyRange with the specified vector<Frequency >.
 	 * @param frequencyRange The vector<Frequency > value to which frequencyRange is to be set.
 	 
 		
 			
 	 */
 	void setFrequencyRange (vector<Frequency > frequencyRange);
  		
	
	
	


	
	// ===> Attribute integrationTime
	
	
	

	
 	/**
 	 * Get integrationTime.
 	 * @return integrationTime as Interval
 	 */
 	Interval getIntegrationTime() const;
	
 
 	
 	
 	/**
 	 * Set integrationTime with the specified Interval.
 	 * @param integrationTime The Interval value to which integrationTime is to be set.
 	 
 		
 			
 	 */
 	void setIntegrationTime (Interval integrationTime);
  		
	
	
	


	
	// ===> Attribute polarizationTypes
	
	
	

	
 	/**
 	 * Get polarizationTypes.
 	 * @return polarizationTypes as vector<PolarizationTypeMod::PolarizationType >
 	 */
 	vector<PolarizationTypeMod::PolarizationType > getPolarizationTypes() const;
	
 
 	
 	
 	/**
 	 * Set polarizationTypes with the specified vector<PolarizationTypeMod::PolarizationType >.
 	 * @param polarizationTypes The vector<PolarizationTypeMod::PolarizationType > value to which polarizationTypes is to be set.
 	 
 		
 			
 	 */
 	void setPolarizationTypes (vector<PolarizationTypeMod::PolarizationType > polarizationTypes);
  		
	
	
	


	
	// ===> Attribute correctionValidity
	
	
	

	
 	/**
 	 * Get correctionValidity.
 	 * @return correctionValidity as bool
 	 */
 	bool getCorrectionValidity() const;
	
 
 	
 	
 	/**
 	 * Set correctionValidity with the specified bool.
 	 * @param correctionValidity The bool value to which correctionValidity is to be set.
 	 
 		
 			
 	 */
 	void setCorrectionValidity (bool correctionValidity);
  		
	
	
	


	
	// ===> Attribute phaseAnt
	
	
	

	
 	/**
 	 * Get phaseAnt.
 	 * @return phaseAnt as vector<float >
 	 */
 	vector<float > getPhaseAnt() const;
	
 
 	
 	
 	/**
 	 * Set phaseAnt with the specified vector<float >.
 	 * @param phaseAnt The vector<float > value to which phaseAnt is to be set.
 	 
 		
 			
 	 */
 	void setPhaseAnt (vector<float > phaseAnt);
  		
	
	
	


	
	// ===> Attribute phaseAntRMS
	
	
	

	
 	/**
 	 * Get phaseAntRMS.
 	 * @return phaseAntRMS as vector<float >
 	 */
 	vector<float > getPhaseAntRMS() const;
	
 
 	
 	
 	/**
 	 * Set phaseAntRMS with the specified vector<float >.
 	 * @param phaseAntRMS The vector<float > value to which phaseAntRMS is to be set.
 	 
 		
 			
 	 */
 	void setPhaseAntRMS (vector<float > phaseAntRMS);
  		
	
	
	


	
	// ===> Attribute amplitudeAnt
	
	
	

	
 	/**
 	 * Get amplitudeAnt.
 	 * @return amplitudeAnt as vector<float >
 	 */
 	vector<float > getAmplitudeAnt() const;
	
 
 	
 	
 	/**
 	 * Set amplitudeAnt with the specified vector<float >.
 	 * @param amplitudeAnt The vector<float > value to which amplitudeAnt is to be set.
 	 
 		
 			
 	 */
 	void setAmplitudeAnt (vector<float > amplitudeAnt);
  		
	
	
	


	
	// ===> Attribute amplitudeAntRMS
	
	
	

	
 	/**
 	 * Get amplitudeAntRMS.
 	 * @return amplitudeAntRMS as vector<float >
 	 */
 	vector<float > getAmplitudeAntRMS() const;
	
 
 	
 	
 	/**
 	 * Set amplitudeAntRMS with the specified vector<float >.
 	 * @param amplitudeAntRMS The vector<float > value to which amplitudeAntRMS is to be set.
 	 
 		
 			
 	 */
 	void setAmplitudeAntRMS (vector<float > amplitudeAntRMS);
  		
	
	
	


	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calDataId
	
	
	

	
 	/**
 	 * Get calDataId.
 	 * @return calDataId as Tag
 	 */
 	Tag getCalDataId() const;
	
 
 	
 	
 	/**
 	 * Set calDataId with the specified Tag.
 	 * @param calDataId The Tag value to which calDataId is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCalDataId (Tag calDataId);
  		
	
	
	


	
	// ===> Attribute calReductionId
	
	
	

	
 	/**
 	 * Get calReductionId.
 	 * @return calReductionId as Tag
 	 */
 	Tag getCalReductionId() const;
	
 
 	
 	
 	/**
 	 * Set calReductionId with the specified Tag.
 	 * @param calReductionId The Tag value to which calReductionId is to be set.
 	 
 		
 			
 	 * @throw IllegalAccessException If an attempt is made to change this field after is has been added to the table.
 	 		
 	 */
 	void setCalReductionId (Tag calReductionId);
  		
	
	
	


	///////////
	// Links //
	///////////
	
	

	
		
	/**
	 * calDataId pointer to the row in the CalData table having CalData.calDataId == calDataId
	 * @return a CalDataRow*
	 * 
	 
	 */
	 CalDataRow* getCalDataUsingCalDataId();
	 

	

	

	
		
	/**
	 * calReductionId pointer to the row in the CalReduction table having CalReduction.calReductionId == calReductionId
	 * @return a CalReductionRow*
	 * 
	 
	 */
	 CalReductionRow* getCalReductionUsingCalReductionId();
	 

	

	
	
	
	/**
	 * Compare each mandatory attribute except the autoincrementable one of this CalAntennaSolutionsRow with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param antennaName
	    
	 * @param atmPhaseCorrection
	    
	 * @param receiverBand
	    
	 * @param basebandName
	    
	 * @param calDataId
	    
	 * @param calReductionId
	    
	 * @param startValidTime
	    
	 * @param endValidTime
	    
	 * @param numReceptor
	    
	 * @param refAntennaName
	    
	 * @param direction
	    
	 * @param frequencyRange
	    
	 * @param integrationTime
	    
	 * @param polarizationTypes
	    
	 * @param correctionValidity
	    
	 * @param phaseAnt
	    
	 * @param phaseAntRMS
	    
	 * @param amplitudeAnt
	    
	 * @param amplitudeAntRMS
	    
	 */ 
	bool compareNoAutoInc(string antennaName, AtmPhaseCorrectionMod::AtmPhaseCorrection atmPhaseCorrection, ReceiverBandMod::ReceiverBand receiverBand, BasebandNameMod::BasebandName basebandName, Tag calDataId, Tag calReductionId, ArrayTime startValidTime, ArrayTime endValidTime, int numReceptor, string refAntennaName, vector<Angle > direction, vector<Frequency > frequencyRange, Interval integrationTime, vector<PolarizationTypeMod::PolarizationType > polarizationTypes, bool correctionValidity, vector<float > phaseAnt, vector<float > phaseAntRMS, vector<float > amplitudeAnt, vector<float > amplitudeAntRMS);
	
	

	
	/**
	 * Compare each mandatory value (i.e. not in the key) attribute  with 
	 * the corresponding parameters and return true if there is a match and false otherwise.
	 	
	 * @param startValidTime
	    
	 * @param endValidTime
	    
	 * @param numReceptor
	    
	 * @param refAntennaName
	    
	 * @param direction
	    
	 * @param frequencyRange
	    
	 * @param integrationTime
	    
	 * @param polarizationTypes
	    
	 * @param correctionValidity
	    
	 * @param phaseAnt
	    
	 * @param phaseAntRMS
	    
	 * @param amplitudeAnt
	    
	 * @param amplitudeAntRMS
	    
	 */ 
	bool compareRequiredValue(ArrayTime startValidTime, ArrayTime endValidTime, int numReceptor, string refAntennaName, vector<Angle > direction, vector<Frequency > frequencyRange, Interval integrationTime, vector<PolarizationTypeMod::PolarizationType > polarizationTypes, bool correctionValidity, vector<float > phaseAnt, vector<float > phaseAntRMS, vector<float > amplitudeAnt, vector<float > amplitudeAntRMS); 
		 
	
	/**
	 * Return true if all required attributes of the value part are equal to their homologues
	 * in x and false otherwise.
	 *
	 * @param x a pointer on the CalAntennaSolutionsRow whose required attributes of the value part 
	 * will be compared with those of this.
	 * @return a boolean.
	 */
	bool equalByRequiredValue(CalAntennaSolutionsRow* x) ;
	
#ifndef WITHOUT_ACS
	/**
	 * Return this row in the form of an IDL struct.
	 * @return The values of this row as a CalAntennaSolutionsRowIDL struct.
	 */
	asdmIDL::CalAntennaSolutionsRowIDL *toIDL() const;
	
	/**
	 * Define the content of a CalAntennaSolutionsRowIDL struct from the values
	 * found in this row.
	 *
	 * @param x a reference to the CalAntennaSolutionsRowIDL struct to be set.
	 *
	 */
	 void toIDL(asdmIDL::CalAntennaSolutionsRowIDL& x) const;
#endif
	
#ifndef WITHOUT_ACS
	/**
	 * Fill the values of this row from the IDL struct CalAntennaSolutionsRowIDL.
	 * @param x The IDL struct containing the values used to fill this row.
	 * @throws ConversionException
	 */
	void setFromIDL (asdmIDL::CalAntennaSolutionsRowIDL x) ;
#endif
	
	/**
	 * Return this row in the form of an XML string.
	 * @return The values of this row as an XML string.
	 */
	std::string toXML() const;

	/**
	 * Fill the values of this row from an XML string 
	 * that was produced by the toXML() method.
	 * @param rowDoc the XML string being used to set the values of this row.
	 * @throws ConversionException
	 */
	void setFromXML (std::string rowDoc) ;

	/// @cond DISPLAY_PRIVATE	
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////

	std::map<std::string, CalAntennaSolutionsAttributeFromBin> fromBinMethods;
void antennaNameFromBin( EndianIStream& eis);
void atmPhaseCorrectionFromBin( EndianIStream& eis);
void receiverBandFromBin( EndianIStream& eis);
void basebandNameFromBin( EndianIStream& eis);
void calDataIdFromBin( EndianIStream& eis);
void calReductionIdFromBin( EndianIStream& eis);
void startValidTimeFromBin( EndianIStream& eis);
void endValidTimeFromBin( EndianIStream& eis);
void numReceptorFromBin( EndianIStream& eis);
void refAntennaNameFromBin( EndianIStream& eis);
void directionFromBin( EndianIStream& eis);
void frequencyRangeFromBin( EndianIStream& eis);
void integrationTimeFromBin( EndianIStream& eis);
void polarizationTypesFromBin( EndianIStream& eis);
void correctionValidityFromBin( EndianIStream& eis);
void phaseAntFromBin( EndianIStream& eis);
void phaseAntRMSFromBin( EndianIStream& eis);
void amplitudeAntFromBin( EndianIStream& eis);
void amplitudeAntRMSFromBin( EndianIStream& eis);

	

	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalAntennaSolutionsTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.
	  */
	 static CalAntennaSolutionsRow* fromBin(EndianIStream& eis, CalAntennaSolutionsTable& table, const std::vector<std::string>& attributesSeq);	 
 
 	 /**
 	  * Parses a string t and assign the result of the parsing to the attribute of name attributeName.
 	  *
 	  * @param attributeName the name of the attribute whose value is going to be defined.
 	  * @param t the string to be parsed into a value given to the attribute of name attributeName.
 	  */
 	 void fromText(const std::string& attributeName, const std::string&  t);
     /// @endcond			

private:
	/**
	 * The table to which this row belongs.
	 */
	CalAntennaSolutionsTable &table;
	/**
	 * Whether this row has been added to the table or not.
	 */
	bool hasBeenAdded;

	// This method is used by the Table class when this row is added to the table.
	void isAdded(bool added);


	/**
	 * Create a CalAntennaSolutionsRow.
	 * <p>
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 */ 
	CalAntennaSolutionsRow (CalAntennaSolutionsTable &table);

	/**
	 * Create a CalAntennaSolutionsRow using a copy constructor mechanism.
	 * <p>
	 * Given a CalAntennaSolutionsRow row and a CalAntennaSolutionsTable table, the method creates a new
	 * CalAntennaSolutionsRow owned by table. Each attribute of the created row is a copy (deep)
	 * of the corresponding attribute of row. The method does not add the created
	 * row to its table, its simply parents it to table, a call to the add method
	 * has to be done in order to get the row added (very likely after having modified
	 * some of its attributes).
	 * If row is null then the method returns a row with default values for its attributes. 
	 *
	 * This constructor is private because only the
	 * table can create rows.  All rows know the table
	 * to which they belong.
	 * @param table The table to which this row belongs.
	 * @param row  The row which is to be copied.
	 */
	 CalAntennaSolutionsRow (CalAntennaSolutionsTable &table, CalAntennaSolutionsRow &row);
	 	
	////////////////////////////////
	// Intrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute antennaName
	
	

	string antennaName;

	
	
 	

	
	// ===> Attribute atmPhaseCorrection
	
	

	AtmPhaseCorrectionMod::AtmPhaseCorrection atmPhaseCorrection;

	
	
 	

	
	// ===> Attribute basebandName
	
	

	BasebandNameMod::BasebandName basebandName;

	
	
 	

	
	// ===> Attribute receiverBand
	
	

	ReceiverBandMod::ReceiverBand receiverBand;

	
	
 	

	
	// ===> Attribute startValidTime
	
	

	ArrayTime startValidTime;

	
	
 	

	
	// ===> Attribute endValidTime
	
	

	ArrayTime endValidTime;

	
	
 	

	
	// ===> Attribute numReceptor
	
	

	int numReceptor;

	
	
 	

	
	// ===> Attribute refAntennaName
	
	

	string refAntennaName;

	
	
 	

	
	// ===> Attribute direction
	
	

	vector<Angle > direction;

	
	
 	

	
	// ===> Attribute frequencyRange
	
	

	vector<Frequency > frequencyRange;

	
	
 	

	
	// ===> Attribute integrationTime
	
	

	Interval integrationTime;

	
	
 	

	
	// ===> Attribute polarizationTypes
	
	

	vector<PolarizationTypeMod::PolarizationType > polarizationTypes;

	
	
 	

	
	// ===> Attribute correctionValidity
	
	

	bool correctionValidity;

	
	
 	

	
	// ===> Attribute phaseAnt
	
	

	vector<float > phaseAnt;

	
	
 	

	
	// ===> Attribute phaseAntRMS
	
	

	vector<float > phaseAntRMS;

	
	
 	

	
	// ===> Attribute amplitudeAnt
	
	

	vector<float > amplitudeAnt;

	
	
 	

	
	// ===> Attribute amplitudeAntRMS
	
	

	vector<float > amplitudeAntRMS;

	
	
 	

	////////////////////////////////
	// Extrinsic Table Attributes //
	////////////////////////////////
	
	
	// ===> Attribute calDataId
	
	

	Tag calDataId;

	
	
 	

	
	// ===> Attribute calReductionId
	
	

	Tag calReductionId;

	
	
 	

	///////////
	// Links //
	///////////
	
	
		

	 

	

	
		

	 

	

	
/*
	////////////////////////////////////////////////////////////
	// binary-deserialization material from an EndianIStream  //
	////////////////////////////////////////////////////////////
	std::map<std::string, CalAntennaSolutionsAttributeFromBin> fromBinMethods;
void antennaNameFromBin( EndianIStream& eis);
void atmPhaseCorrectionFromBin( EndianIStream& eis);
void receiverBandFromBin( EndianIStream& eis);
void basebandNameFromBin( EndianIStream& eis);
void calDataIdFromBin( EndianIStream& eis);
void calReductionIdFromBin( EndianIStream& eis);
void startValidTimeFromBin( EndianIStream& eis);
void endValidTimeFromBin( EndianIStream& eis);
void numReceptorFromBin( EndianIStream& eis);
void refAntennaNameFromBin( EndianIStream& eis);
void directionFromBin( EndianIStream& eis);
void frequencyRangeFromBin( EndianIStream& eis);
void integrationTimeFromBin( EndianIStream& eis);
void polarizationTypesFromBin( EndianIStream& eis);
void correctionValidityFromBin( EndianIStream& eis);
void phaseAntFromBin( EndianIStream& eis);
void phaseAntRMSFromBin( EndianIStream& eis);
void amplitudeAntFromBin( EndianIStream& eis);
void amplitudeAntRMSFromBin( EndianIStream& eis);

	
*/
	
	///////////////////////////////////
	// text-deserialization material //
	///////////////////////////////////
	std::map<std::string, CalAntennaSolutionsAttributeFromText> fromTextMethods;
	
void antennaNameFromText (const string & s);
	
	
void atmPhaseCorrectionFromText (const string & s);
	
	
void receiverBandFromText (const string & s);
	
	
void basebandNameFromText (const string & s);
	
	
void calDataIdFromText (const string & s);
	
	
void calReductionIdFromText (const string & s);
	
	
void startValidTimeFromText (const string & s);
	
	
void endValidTimeFromText (const string & s);
	
	
void numReceptorFromText (const string & s);
	
	
void refAntennaNameFromText (const string & s);
	
	
void directionFromText (const string & s);
	
	
void frequencyRangeFromText (const string & s);
	
	
void integrationTimeFromText (const string & s);
	
	
void polarizationTypesFromText (const string & s);
	
	
void correctionValidityFromText (const string & s);
	
	
void phaseAntFromText (const string & s);
	
	
void phaseAntRMSFromText (const string & s);
	
	
void amplitudeAntFromText (const string & s);
	
	
void amplitudeAntRMSFromText (const string & s);
	

		
	
	/**
	 * Serialize this into a stream of bytes written to an EndianOSStream.
	 * @param eoss the EndianOSStream to be written to
	 */
	 void toBin(EndianOSStream& eoss);
	 	 
	 /**
	  * Deserialize a stream of bytes read from an EndianIStream to build a PointingRow.
	  * @param eiss the EndianIStream to be read.
	  * @param table the CalAntennaSolutionsTable to which the row built by deserialization will be parented.
	  * @param attributesSeq a vector containing the names of the attributes . The elements order defines the order 
	  * in which the attributes are written in the binary serialization.

	 static CalAntennaSolutionsRow* fromBin(EndianIStream& eis, CalAntennaSolutionsTable& table, const std::vector<std::string>& attributesSeq);	 
		*/
};

} // End namespace asdm

#endif /* CalAntennaSolutions_CLASS */
