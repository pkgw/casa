//# tFlagAgentDisplay.cc This file contains the unit tests of the FlagAgentBase class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#include <flagging/Flagging/FlagAgentDisplay.h>
#include <flagging/Flagging/FlagAgentManual.h>
#include <iostream>

using namespace casacore;
using namespace casa;

void deleteFlags(string inputFile,Record dataSelection,vector<Record> agentParameters)
{
	// Some test execution info
	cout << "STEP 1: CLEAN FLAGS ..." << endl;

	// Stats variables declaration
	unsigned long nBuffers = 0;
	unsigned long cumRows = 0;
	timeval start,stop;
	double elapsedTime = 0;
	Double ntime = 0;

	// Start clock
	gettimeofday(&start,0);

	// Parse time interval
	if (dataSelection.fieldNumber ("ntime") >= 0)
	{
		ntime = dataSelection.asDouble("ntime");
	}

	// Create Flag Data Handler
	FlagDataHandler *dh = new FlagMSHandler(inputFile,FlagDataHandler::COMPLETE_SCAN_UNMAPPED,ntime);

	// Enable profiling in the Flag Data Handler
	dh->setProfiling(false);

	// Open MS
	dh->open();

	// Parse data selection to Flag Data Handler
	dh->setDataSelection(dataSelection);

	// Select data (creating selected MS)
	dh->selectData();

	// Generate iterators and vis buffers
	dh->generateIterator();

	// Create agent list
	Int agentNumber = 1;
	FlagAgentList agentList;
	FlagAgentManual *flaggingAgent = NULL;
	for (vector<Record>::iterator iter=agentParameters.begin();iter != agentParameters.end();iter++)
	{
		stringstream agentName;
		agentName << agentNumber;
		iter->define("name","FlagAgentManual_" + agentName.str());
		flaggingAgent = new FlagAgentManual(dh,*iter,false,false);
		agentList.push_back(flaggingAgent);
		agentNumber++;
	}

	// Enable profiling in the Flag Agent
	agentList.setProfiling(false);

	// Start Flag Agent
	agentList.start();

	// Set cout precision
	cout.precision(20);

	// iterate over chunks
	while (dh->nextChunk())
	{
		// iterate over visBuffers
		while (dh->nextBuffer())
		{
			//cout << "Chunk:" << dh->chunkNo << " " << "Buffer:" << dh->bufferNo << " ";
			nBuffers += 1;
/*
			if (dh->visibilityBuffer_p->get()->observationId().nelements() > 1)
			{
				cout << "Observation:"
					 << dh->visibilityBuffer_p->get()->observationId()[0] << "~"
					 << dh->visibilityBuffer_p->get()->observationId()[dh->visibilityBuffer_p->get()->observationId().nelements()-1] << " ";
			}
			else
			{
				cout << "Observation:" << dh->visibilityBuffer_p->get()->observationId()[0] << " ";
			}

			cout << "Array:" << dh->visibilityBuffer_p->get()->arrayId() << " ";

			if (dh->visibilityBuffer_p->get()->scan().nelements() > 1)
			{
				cout << "Scan:"
					 << dh->visibilityBuffer_p->get()->scan()[0] << "~"
					 << dh->visibilityBuffer_p->get()->scan()[dh->visibilityBuffer_p->get()->scan().nelements()-1] << " ";
			}
			else
			{
				cout << "Scan:" << dh->visibilityBuffer_p->get()->scan()[0] << " ";
			}

			cout << "Field:" << dh->visibilityBuffer_p->get()->fieldId() << " " ;

			cout << "Spw:" << dh->visibilityBuffer_p->get()->spectralWindow() << " ";

			if (dh->visibilityBuffer_p->get()->time().nelements() > 1)
			{
				cout << "Time:"
					 << dh->visibilityBuffer_p->get()->time()[0] << "~"
					 << dh->visibilityBuffer_p->get()->time()[dh->visibilityBuffer_p->get()->time().nelements()-1] << " ";
			}
			else
			{
				cout << "Time:" << dh->visibilityBuffer_p->get()->time()[0] << " ";
			}

			if (dh->visibilityBuffer_p->get()->antenna1().nelements() > 1)
			{
				cout << "Antenna1:"
					 << dh->visibilityBuffer_p->get()->antenna1()[0] << "~"
					 << dh->visibilityBuffer_p->get()->antenna1()[dh->visibilityBuffer_p->get()->antenna1().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna1:" << dh->visibilityBuffer_p->get()->antenna1()[0] << " ";
			}

			if (dh->visibilityBuffer_p->get()->antenna2().nelements() > 1)
			{
				cout << "Antenna2:"
					 << dh->visibilityBuffer_p->get()->antenna2()[0] << "~"
					 << dh->visibilityBuffer_p->get()->antenna2()[dh->visibilityBuffer_p->get()->antenna2().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna2:" << dh->visibilityBuffer_p->get()->antenna2()[0] << " ";
			}
*/
			//cout << "nRows:" << dh->visibilityBuffer_p->get()->nRow() <<endl;
			cumRows += dh->visibilityBuffer_p->get()->nRow();

			// Apply flags
			agentList.apply();

			// Flush flags to MS
			dh->flushFlags();
		}

		// Print stats from each agent
		agentList.chunkSummary();
	}

	// Print total stats from each agent
	agentList.tableSummary();

	// Stop Flag Agent
	agentList.terminate();
	agentList.join();

	// Close MS
	dh->close();

	// Clear Flag Agent List
	agentList.clear();

	// Delete Flag Data Handler (delete VisBuffer, therefore stop VLAT)
	delete dh;

	// Stop clock
	gettimeofday(&stop,0);
	elapsedTime = (stop.tv_sec-start.tv_sec)*1000.0+(stop.tv_usec-start.tv_usec)/1000.0;

	// Report elapsed time
	cout << "Total Time [s]:" << elapsedTime/1000.0 << " Total number of rows:" << cumRows <<" Total number of Buffers:" << nBuffers <<endl;

}

void writeFlags(string inputFile,Record dataSelection,vector<Record> agentParameters)
{
	// Some test execution info
	cout << "STEP 2: WRITE FLAGS ..." << endl;

	// Stats variables declaration
	unsigned long nBuffers = 0;
	unsigned long cumRows = 0;
	timeval start,stop;
	double elapsedTime = 0;
	Double ntime = 0;

	// Start clock
	gettimeofday(&start,0);

	// Parse time interval
	if (dataSelection.fieldNumber ("ntime") >= 0)
	{
		ntime = dataSelection.asDouble("ntime");
	}

	// Create Flag Data Handler
	FlagDataHandler *dh = new FlagMSHandler(inputFile,FlagDataHandler::COMPLETE_SCAN_MAP_ANTENNA_PAIRS_ONLY,ntime);

	// Enable profiling in the Flag Data Handler
	dh->setProfiling(false);

	// Open MS
	dh->open();

	// Set data selection
	dh->setDataSelection(dataSelection);

	// Select data (creating selected MS)
	dh->selectData();

	// Generate iterators and vis buffers
	dh->generateIterator();

	// Create agent list :  One ManualFlag, and One Display
	Int agentNumber = 1;
	FlagAgentList agentList;
	FlagAgentDisplay *flaggingAgent = NULL;
        FlagAgentManual *flaggingAgent2 = NULL;
	//        FlagAgentSummary *flaggingAgent3 = NULL;

                
        // Add a manual-flag agent (the list has only one here)
	for (vector<Record>::iterator iter=agentParameters.begin();iter != agentParameters.end();iter++)
	{
          // (1) Change the spw parameter for manual flag only. Not for the display
          String spwstr("");
          if( iter->fieldNumber("spw") >=0 )  
            {
              spwstr = iter->asString("spw");
              cout << "Original SPW string : " << spwstr << endl;
            }
          iter->define("spw","*:35~40");
          // (2) Make the manual-flag agent 
          stringstream agentName;
          agentName << agentNumber;
          iter->define("name","FlagAgentManual_" + agentName.str());
          flaggingAgent2 = new FlagAgentManual(dh,*iter,false,true);
          agentList.push_back(flaggingAgent2);
          agentNumber++;
          // (3) Change the spw parameter back to what it was...
          if( spwstr == "" ) iter->removeField( iter->fieldNumber("spw") );
          else iter->define("spw",spwstr);
	}
        cout << "Agent Count after Manual flag only : " << agentNumber-1 << endl;

	/*	
	// Add a summary agent
	for (vector<Record>::iterator iter=agentParameters.begin();iter != agentParameters.end();iter++)
	{
		stringstream agentName;
		agentName << agentNumber;
		iter->define("name","FlagAgentSummary_" + agentName.str());
		flaggingAgent3 = new FlagAgentSummary(dh,*iter);
		agentList.push_back(flaggingAgent3);
		agentNumber++;
	}
        cout << "Agent Count after adding Summary : " << agentNumber-1 << endl;
	*/

        // Add the display agent
	for (vector<Record>::iterator iter=agentParameters.begin();iter != agentParameters.end();iter++)
	{
		stringstream agentName;
		agentName << agentNumber;
		iter->define("name","FlagAgentDisplay_" + agentName.str());
		flaggingAgent = new FlagAgentDisplay(dh,*iter,false/*writePrivFlagCube*/);
		agentList.push_back(flaggingAgent);
		agentNumber++;
	}
        cout << "Agent Count after Display too : " << agentNumber-1 << endl;

	// Enable profiling mode
	agentList.setProfiling(false);

	// Start Flag Agent
	agentList.start();

	// Set cout precision
	cout.precision(20);

	// iterate over chunks
	while (dh->nextChunk())
	{
		// iterate over visBuffers
		while (dh->nextBuffer())
		{
			//cout << "Chunk:" << dh->chunkNo << " " << "Buffer:" << dh->bufferNo << " ";
			nBuffers += 1;
/*
			if (dh->visibilityBuffer_p->get()->observationId().nelements() > 1)
			{
				cout << "Observation:"
					 << dh->visibilityBuffer_p->get()->observationId()[0] << "~"
					 << dh->visibilityBuffer_p->get()->observationId()[dh->visibilityBuffer_p->get()->observationId().nelements()-1] << " ";
			}
			else
			{
				cout << "Observation:" << dh->visibilityBuffer_p->get()->observationId()[0] << " ";
			}

			cout << "Array:" << dh->visibilityBuffer_p->get()->arrayId() << " ";

			if (dh->visibilityBuffer_p->get()->scan().nelements() > 1)
			{
				cout << "Scan:"
					 << dh->visibilityBuffer_p->get()->scan()[0] << "~"
					 << dh->visibilityBuffer_p->get()->scan()[dh->visibilityBuffer_p->get()->scan().nelements()-1] << " ";
			}
			else
			{
				cout << "Scan:" << dh->visibilityBuffer_p->get()->scan()[0] << " ";
			}

			cout << "Field:" << dh->visibilityBuffer_p->get()->fieldId() << " " ;

			cout << "Spw:" << dh->visibilityBuffer_p->get()->spectralWindow() << " ";

			if (dh->visibilityBuffer_p->get()->time().nelements() > 1)
			{
				cout << "Time:"
					 << dh->visibilityBuffer_p->get()->time()[0] << "~"
					 << dh->visibilityBuffer_p->get()->time()[dh->visibilityBuffer_p->get()->time().nelements()-1] << " ";
			}
			else
			{
				cout << "Time:" << dh->visibilityBuffer_p->get()->time()[0] << " ";
			}

			if (dh->visibilityBuffer_p->get()->antenna1().nelements() > 1)
			{
				cout << "Antenna1:"
					 << dh->visibilityBuffer_p->get()->antenna1()[0] << "~"
					 << dh->visibilityBuffer_p->get()->antenna1()[dh->visibilityBuffer_p->get()->antenna1().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna1:" << dh->visibilityBuffer_p->get()->antenna1()[0] << " ";
			}

			if (dh->visibilityBuffer_p->get()->antenna2().nelements() > 1)
			{
				cout << "Antenna2:"
					 << dh->visibilityBuffer_p->get()->antenna2()[0] << "~"
					 << dh->visibilityBuffer_p->get()->antenna2()[dh->visibilityBuffer_p->get()->antenna2().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna2:" << dh->visibilityBuffer_p->get()->antenna2()[0] << " ";
			}
*/
			//cout << "nRows:" << dh->visibilityBuffer_p->get()->nRow() <<endl;
			cumRows += dh->visibilityBuffer_p->get()->nRow();

			// IMPORTANT: Pre-load vis cube, antenna1 and antenna2
			dh->visibilityBuffer_p->get()->visCube();

			// Apply flags
			agentList.apply();

			// Flush flags to MS
			dh->flushFlags();
		}

		// Print stats from each agent
		agentList.chunkSummary();
	}

	// Print total stats from each agent
	agentList.tableSummary();

	// Stop Flag Agent
	agentList.terminate();
	agentList.join();

        // Gather reports from all Agents that make them
        FlagReport combinedReport = agentList.gatherReports();

	// Print the combined Record (for debugging)
	stringstream replist;
	combinedReport.print(replist);
        cout << " Combined Report : " << endl << replist.str() << endl;

	// Send them to display.
        flaggingAgent->displayReports(combinedReport);


	// Close MS
	dh->close();

	// Clear Flag Agent List
	agentList.clear();

	// Delete Flag Data Handler (delete VisBuffer, therefore stop VLAT)
	delete dh;

	// Stop clock
	gettimeofday(&stop,0);
	elapsedTime = (stop.tv_sec-start.tv_sec)*1000.0+(stop.tv_usec-start.tv_usec)/1000.0;

	// Report elapsed time
	cout << "Total Time [s]:" << elapsedTime/1000.0 << " Total number of rows:" << cumRows <<" Total number of Buffers:" << nBuffers <<endl;

}

bool checkFlags(string targetFile,string referenceFile, Record dataSelection)
{
	// Some test execution info
	cout << "STEP 3: CHECK FLAGS ..." << endl;

	// Declare variables
	Cube<Bool> targetFileFlags;
	Cube<Bool> referenceFileFlags;
	IPosition targetFileFlagsShape;
	IPosition referenceFileFlagsShape;

	// Execution control variables declaration
	bool returnCode=true;

	// Create data handler for target file
	FlagDataHandler *targetFiledh = new FlagMSHandler(targetFile,FlagDataHandler::COMPLETE_SCAN_UNMAPPED);
	targetFiledh->open();
	targetFiledh->setDataSelection(dataSelection);
	targetFiledh->selectData();
	targetFiledh->generateIterator();

	// Create data handler for reference file
	FlagDataHandler *referenceFiledh = new FlagMSHandler(referenceFile,FlagDataHandler::COMPLETE_SCAN_UNMAPPED);
	referenceFiledh->open();
	referenceFiledh->setDataSelection(dataSelection);
	referenceFiledh->selectData();
	referenceFiledh->generateIterator();

	// iterate over chunks
	while (targetFiledh->nextChunk() and referenceFiledh->nextChunk())
	{
		// iterate over visBuffers
		while (targetFiledh->nextBuffer() and referenceFiledh->nextBuffer())
		{
			cout << "Chunk:" << targetFiledh->chunkNo << " " << "Buffer:" << targetFiledh->bufferNo << " ";

			if (targetFiledh->visibilityBuffer_p->get()->observationId().nelements() > 1)
			{
				cout << "Observation:"
					 << targetFiledh->visibilityBuffer_p->get()->observationId()[0] << "~"
					 << targetFiledh->visibilityBuffer_p->get()->observationId()[targetFiledh->visibilityBuffer_p->get()->observationId().nelements()-1] << " ";
			}
			else
			{
				cout << "Observation:" << targetFiledh->visibilityBuffer_p->get()->observationId()[0] << " ";
			}

			cout << "Array:" << targetFiledh->visibilityBuffer_p->get()->arrayId() << " ";

			if (targetFiledh->visibilityBuffer_p->get()->scan().nelements() > 1)
			{
				cout << "Scan:"
					 << targetFiledh->visibilityBuffer_p->get()->scan()[0] << "~"
					 << targetFiledh->visibilityBuffer_p->get()->scan()[targetFiledh->visibilityBuffer_p->get()->scan().nelements()-1] << " ";
			}
			else
			{
				cout << "Scan:" << targetFiledh->visibilityBuffer_p->get()->scan()[0] << " ";
			}

			cout << "Field:" << targetFiledh->visibilityBuffer_p->get()->fieldId() << " " ;

			cout << "Spw:" << targetFiledh->visibilityBuffer_p->get()->spectralWindow() << " ";

			if (targetFiledh->visibilityBuffer_p->get()->time().nelements() > 1)
			{
				cout << "Time:"
					 << targetFiledh->visibilityBuffer_p->get()->time()[0] << "~"
					 << targetFiledh->visibilityBuffer_p->get()->time()[targetFiledh->visibilityBuffer_p->get()->time().nelements()-1] << " ";
			}
			else
			{
				cout << "Time:" << targetFiledh->visibilityBuffer_p->get()->time()[0] << " ";
			}

			if (targetFiledh->visibilityBuffer_p->get()->antenna1().nelements() > 1)
			{
				cout << "Antenna1:"
					 << targetFiledh->visibilityBuffer_p->get()->antenna1()[0] << "~"
					 << targetFiledh->visibilityBuffer_p->get()->antenna1()[targetFiledh->visibilityBuffer_p->get()->antenna1().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna1:" << targetFiledh->visibilityBuffer_p->get()->antenna1()[0] << " ";
			}

			if (targetFiledh->visibilityBuffer_p->get()->antenna2().nelements() > 1)
			{
				cout << "Antenna2:"
					 << targetFiledh->visibilityBuffer_p->get()->antenna2()[0] << "~"
					 << targetFiledh->visibilityBuffer_p->get()->antenna2()[targetFiledh->visibilityBuffer_p->get()->antenna2().nelements()-1] << " ";
			}
			else
			{
				cout << "Antenna2:" << targetFiledh->visibilityBuffer_p->get()->antenna2()[0] << " ";
			}

			cout << "nRows:" << targetFiledh->visibilityBuffer_p->get()->nRow() <<endl;

			targetFileFlags = targetFiledh->visibilityBuffer_p->get()->flagCube();
			referenceFileFlags = referenceFiledh->visibilityBuffer_p->get()->flagCube();
			IPosition targetFileFlagsShape = targetFileFlags.shape();
			IPosition referenceFileFlagsShape = referenceFileFlags.shape();

			if (targetFileFlagsShape != referenceFileFlagsShape)
			{
				cerr << "Target and reference flag cubes have different shape " << endl;
				returnCode = false;
			}

			for (Int row=0;row<targetFileFlagsShape(3);row++)
			{
				for (Int chan=0;chan<targetFileFlagsShape(2);chan++)
				{
					for (Int corr=0;corr<targetFileFlagsShape(1);corr++)
					{
						if (targetFileFlags(corr,chan,row) != referenceFileFlags(corr,chan,row))
						{
							cerr << "Flags for Chunk=" << targetFiledh->chunkNo
									<< " buffer=" << targetFiledh->bufferNo
									<< " row=" << row
									<< " chan=" << chan
									<< " corr=" << corr
									<< ", are incorrect, calculated flags=" << targetFileFlags(corr,chan,row)
									<< " reference flags=" << referenceFileFlags(corr,chan,row) << endl;
							returnCode = false;
						}
					}
				}
			}
		}
	}

	delete targetFiledh;
	delete referenceFiledh;

	return returnCode;
}

int main(int argc, char **argv)
{
	// Parsing variables declaration
	string parameter, value;
	string targetFile,referenceFile;
	string array,scan,timerange,field,spw,antenna,uvrange,correlation,observation,intent;
	string expression,datacolumn,ntime;

	// Execution control variables declaration
	bool deleteFlagsActivated=false;
	bool checkFlagsActivated=false;
	bool returnCode=true;

	// Parse input parameters
	Record agentParameters;
	Record dataSelection;
        vector<Record> recordlist;
	for (unsigned short i=0;i<argc-1;i++)
	{
		parameter = string(argv[i]);
		value = string(argv[i+1]);

		if (parameter == string("-targetFile"))
		{
			targetFile = value;
			cout << "Target file is: " << targetFile << endl;
		}
		else if (parameter == string("-referenceFile"))
		{
			referenceFile = value;
			checkFlagsActivated = true;
			cout << "Reference file is: " << referenceFile << endl;
		}
		else if (parameter == string("-unflag"))
		{
			if (value.compare("true") == 0)
			{
				deleteFlagsActivated = true;
				cout << "Clean flags step activated" << endl;
			}
		}
		else if (parameter == string("-array"))
		{
			array = casa::String(value);
			dataSelection.define ("array", array);
			cout << "Array selection is: " << array << endl;
		}
		else if (parameter == string("-scan"))
		{
			scan = casa::String(value);
			dataSelection.define ("scan", scan);
			cout << "Scan selection is: " << scan << endl;
		}
		else if (parameter == string("-timerange"))
		{
			timerange = casa::String(value);
			dataSelection.define ("timerange", timerange);
			cout << "Time range selection is: " << timerange << endl;
		}
		else if (parameter == string("-field"))
		{
			field = casa::String(value);
			dataSelection.define ("field", field);
			cout << "Field selection is: " << field << endl;
		}
		else if (parameter == string("-spw"))
		{
			spw = casa::String(value);
			dataSelection.define ("spw", spw);
			cout << "SPW selection is: " << spw << endl;
		}
		else if (parameter == string("-antenna"))
		{
			antenna = casa::String(value);
			dataSelection.define("antenna",antenna);
			cout << "Antenna selection is: " << antenna << endl;
		}
		else if (parameter == string("-uvrange"))
		{
			uvrange = casa::String(value);
			dataSelection.define ("uvrange", uvrange);
			cout << "UV range selection is: " << uvrange << endl;
		}
		else if (parameter == string("-observation"))
		{
			observation = casa::String(value);
			dataSelection.define ("observation", observation);
			cout << "Observation selection is: " << observation << endl;
		}
		else if (parameter == string("-intent"))
		{
			intent = casa::String(value);
			dataSelection.define ("intent", intent);
			cout << "Scan intention selection is: " << intent << endl;
		}
		else if (parameter == string("-ntime"))
		{
			ntime = casa::String(value);
			dataSelection.define ("ntime", atof(ntime.c_str()));
			cout << "ntime is: " << ntime << endl;
		}
		else if (parameter == string("-correlation"))
		{
			correlation = casa::String(value);
			agentParameters.define ("correlation", correlation);
			cout << "Correlation range selection is: " << correlation << endl;
		}
		else if (parameter == string("-expression"))
		{
			expression = casa::String(value);
			agentParameters.define ("expression", expression);
			cout << "expression is: " << expression << endl;
		}
		else if (parameter == string("-datacolumn"))
		{
			datacolumn = casa::String(value);
			agentParameters.define ("datacolumn", datacolumn);
			cout << "datacolumn is: " << datacolumn << endl;
		}
		else if (parameter == string("-pause"))
		{
                        Bool pause=true;
                        if(value.compare("false")==0) pause=false;
			agentParameters.define ("pause", pause);
			cout << "pause is: " << pause << endl;
		}
		else if (parameter == string("-datadisplay"))
		{
                        Bool datadisplay=true;
                        if(value.compare("false")==0) datadisplay=false;
			agentParameters.define ("datadisplay", datadisplay);
			cout << "datadisplay is: " << datadisplay << endl;
		}
		else if (parameter == string("-reportdisplay"))
		{
                        Bool reportdisplay=true;
                        if(value.compare("false")==0) reportdisplay=false;
			agentParameters.define ("reportdisplay", reportdisplay);
			cout << "reportdisplay is: " << reportdisplay << endl;
		}
		else if (parameter == string("-format"))
		{
		        String format = casa::String(value);
			agentParameters.define ("format", format);
			cout << "format is: " << format << endl;
		}
		else if (parameter == string("-flagspw"))
		{
			spw = casa::String(value);
			dataSelection.define ("spw", spw);
			cout << "SPW selection is: " << spw << endl;
		}

	}


	Record agentParameters_i;
	vector<Record> agentParamersList;

	{
		agentParamersList.push_back(agentParameters);
	}

	if (deleteFlagsActivated) deleteFlags(targetFile,dataSelection,agentParamersList);
	writeFlags(targetFile,dataSelection,agentParamersList);
	if (checkFlagsActivated) returnCode = checkFlags(targetFile,referenceFile,dataSelection);

	if (returnCode)
	{
		exit(0);
	}
	else
	{
		exit(-1);
	}
}
