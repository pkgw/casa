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
 * File NumberFormatException.h
 */

#ifndef NumberFormatException_CLASS
#define NumberFormatException_CLASS

#include <alma/ASDM/InvalidArgumentException.h>

#include <string>

namespace asdm {

/**
 * The NumberFormatException class represents an exception when 
 * an error occurs in converting a numeric value from a string.
 */
  class NumberFormatException : public InvalidArgumentException {
    
  public:
    /**
     * An empty CTOR
     */
    NumberFormatException();

    /**
     * A CTOR with a message.
     */
    NumberFormatException(const std::string &s);

    /**
     * the DTOR
     */
    virtual ~NumberFormatException();
    
  };
} // End namespace asdm

#endif /* NumberFormatException_CLASS */
