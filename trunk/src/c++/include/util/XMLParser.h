#ifndef UTIL_XML_PARSER_H
#define UTIL_XML_PARSER_H

#include <string>
#include <sstream>
#include <vector>
#include "error/internal_error.h"
#include "util/StringUtils.h"
#include "libxml++/libxml++.h"

namespace SAGE {
namespace UTIL {

/// \brief Represents the name of the node as it lies within the hierarchy of an XML document.
///
class QualifiedNodeName
{
public:
  ///
  /// Constructor.
  /// \param name_list The fully qualified path of the node name, comma separated (ie:
  /// "foo,bar,thing".
  /// \throws std::exception if the name_list is empty
  explicit QualifiedNodeName(const std::string & name_list);

  ///
  /// Copy constructor.
  /// \param other The object to copy
  QualifiedNodeName(const QualifiedNodeName & other);

  ///
  /// Renders the stored fully qualified node name as a comma-separated list
  std::string toString() const;

  ///
  /// Returns the number of node names comprising the fully qualified node name this
  /// object represents.
  int getNodeCount() const;

  ///
  /// Returns the first element (the very beginning node name) of the fully qualified
  /// node name.
  const std::string & getFirstName() const;

  ///
  /// Returns the last element (the terminal node name) of the fully qualified
  /// node name.
  const std::string & getLastName() const;

  ///
  /// Returns a copy of this QualifiedNodeName, without the first element.
  /// \throws std::exception if getNodeCount() == 1
  QualifiedNodeName popFront() const;


private:
  QualifiedNodeName& operator= (const QualifiedNodeName & other);

  std::vector<std::string> my_names;
};

inline
QualifiedNodeName::QualifiedNodeName(const std::string & name_list)
{
  // Split up the node name:        
  StringUtils::splitString(name_list, ",", my_names);

  if(my_names.size() == 0)
    throw std::exception();
}

inline
QualifiedNodeName::QualifiedNodeName(const QualifiedNodeName & other)
{
  my_names = other.my_names;
}

inline std::string 
QualifiedNodeName::toString() const
{
  std::ostringstream s;
  
  for(size_t i = 0; i < my_names.size(); ++i)
    s << (i ? "," : "") << my_names[i];
    
  return s.str();
}

inline int 
QualifiedNodeName::getNodeCount() const
{
  return my_names.size();
}

inline const std::string & 
QualifiedNodeName::getFirstName() const
{
  return my_names[0];
}

inline const std::string & 
QualifiedNodeName::getLastName() const
{
  return my_names[my_names.size() - 1];
}

inline QualifiedNodeName 
QualifiedNodeName::popFront() const
{
  if(getNodeCount() < 2)
    throw std::exception();
    
  std::ostringstream s;
      
  for(size_t i = 1; i < my_names.size(); ++i)
    s << (i > 1 ? "," : "") << my_names[i];
                
  return QualifiedNodeName(s.str());
}

/// \brief Provides a helpful wrapper interface to xmlpp for a single XML document
///
/// \par Introduction
///
/// libxml++ is somewhat awkward to use. It works well, but it uses a lot of pointers
/// and so on. Consequently you run into weird memory allocation issues.
///
/// The purpose of the XMLParser class is to provide a wrapper interface for libxml++'s
/// abilities. You can 
class XMLParser
{
public:

  /// A const Node pointer (see xmlpp for reference info)
  typedef const xmlpp::Node * NodePtr;

  /// A vector of xmlpp::Node pointers.
  typedef std::vector<NodePtr> NodePtrVector;

  /// The attribute list.
  typedef xmlpp::Element::AttributeList AttributeList;

  /// @name Well-formed verification
  //@{
  
    ///
    /// Verifies that the given filename is available and well-formed (as an XML document).
    static bool isFileWellFormed(const std::string & filename);
  
    ///
    /// Verifies that the given string is available and well-formed (as an XML document).
    static bool isStringWellFormed(const std::string & xml_text);

  //@}

  /// @name Constructor
  //@{
  
    ///
    /// Constructor.
    /// \param s The name of the XML file or XML text to read in
    /// \param treat_as_file If \c true, \c s is interpreted as a filename; if \c false, \c s
    /// is interpreted as the XML content string itself.
    /// \throws A std::exception if the filename is malformed; you should use
    /// isFileWellFormed() before constructing this object.
    XMLParser(const std::string & s, bool treat_as_file);
  
  //@}

  /// @name Fetching & processing nodes
  //@{
  
    ///
    /// This pulls out the actual text content of an xmlpp::Node.
    /// \param node The node from which text content will be extracted
    /// \returns The text content, or a blank string if there is none
    std::string extractTextContent(NodePtr node) const;
      
    ///
    /// Searches the given XML file for a specific node with text content.
    /// \param node_name The fully qualified name of the node, as in "Foo:Bar:Thing:path"
    /// \param node The node in which to search; if not specified (NULL), starts with the root node of
    /// this object.
    /// \returns \c true if the named node exists and has text content, \c false otherwise
    bool hasNode(const QualifiedNodeName & qnodename, NodePtr node = NULL) const;

    ///
    /// Finds the node.
    ///
    /// Please note: If the node if found, returns the first element found matching the node name. Consider 
    /// the following XML file:
    ///
    /// \verbatim
    /// <Foo>
    ///   <Bar/>
    ///   <Bar>Thing</Bar>
    /// </Foo>
    /// \endverbatim
    ///
    /// A call to getNode(QualifiedNodeName("Foo,Bar")) will return the \b first 'Bar' element, not the
    /// second. Consequently, extractTextContent(getNode(QualifiedNodeName("Foo,Bar"))) will return
    /// an \b empty string, not 'Thing'. You should generally only use this function if you only expect
    /// there to be one instance of the node you're looking for.
    /// 
    /// \param qnodename The fully qualified node name to search for
    /// \param node The node in which to search; if not specified (NULL), starts with the root node of
    /// this object.
    /// \returns The found node
    /// \throws std::exception if the node is not found
    NodePtr getNode(const QualifiedNodeName & qnodename, NodePtr node = NULL) const;

    ///
    /// Finds the all instances of the given node name.
    ///
    /// Consider the following XML file:
    ///
    /// \verbatim
    /// <Foo>
    ///   <Bar/>
    ///   <Bar>Thing</Bar>
    /// </Foo>
    /// \endverbatim
    ///
    /// A call to getNodes(QualifiedNodeName("Foo,Bar")) will return a vector with two elements (both 'Bar' elements).
    ///
    /// \param qnodename The fully qualified node name to search for
    /// \param node The node in which to search; if not specified (NULL), starts with the root node of
    /// this object.
    /// \returns The found node
    /// \throws std::exception if the node is not found
    NodePtrVector getNodes(const QualifiedNodeName & qnodename, NodePtr node = NULL) const;

    ///
    /// Returns the underlying root node from libxml++.
    NodePtr getRootNode() const;

  //@}

  /// @name Fetching & processing node attributes
  //@{
  
    ///
    /// Returns the AttributeList for the given node.
    /// \param node The node in question
    AttributeList getAttributes(NodePtr node) const;

  //@}

  /// @name Outputting contents
  //@{

    ///
    /// Returns the contents of this XML object as an XML string.
    std::string toString() const;

  //@}

private:

  /// @name Disallowed functions:
  //@{
  
    XMLParser(const XMLParser &);
    XMLParser& operator= (const XMLParser &);

  //@}

    xmlpp::DomParser my_parser;
};

//=================================================
//
//  isFileWellFormed(...)
//
//=================================================
inline bool 
XMLParser::isFileWellFormed(const std::string & filename)
{
  try
  {
    XMLParser x(filename, true);
    
    return true;
  }
  catch(const std::exception & e)
  {
    return false;
  }
}

//=================================================
//
//  isStringWellFormed(...)
//
//=================================================
inline bool 
XMLParser::isStringWellFormed(const std::string & xml_text)
{
  try
  {
    XMLParser x(xml_text, false);
    
    return true;
  }
  catch(const std::exception & e)
  {
    return false;
  }
}
            
//====================================================
//
//  CONSTRUCTOR
//
//====================================================
inline
XMLParser::XMLParser(const std::string & s, bool treat_as_filename)
{
  my_parser.set_substitute_entities();

  if(treat_as_filename)
    my_parser.parse_file(s);
  else
    my_parser.parse_memory(s);

  if(getRootNode() == NULL)
    throw std::exception();
}

//====================================================
//
//  hasNode(...)
//
//====================================================
inline bool 
XMLParser::hasNode(const QualifiedNodeName & qnodename, NodePtr node) const
{
  try
  {  
    getNode(qnodename, node);
    
    return true;
  }
  catch(const std::exception & e)
  {
    return false;
  }                                    
}
    
//==========================================================
//
//  getNode(...)
//
//==========================================================
inline XMLParser::NodePtr 
XMLParser::getNode(const QualifiedNodeName & qnodename, NodePtr _node) const
{
  // If the default argument (NULL) was used for _node, then use this object's root node:
  
  NodePtr node = _node == NULL ? getRootNode() : _node;

  // Is the first name correct?

  if(qnodename.getNodeCount() == 1)
  {
    if(node->get_name() == qnodename.getFirstName())
    {
      return node;
    }
    else
    {
      throw std::exception();
    }
  }
  else // nodecount > 1
  {
    QualifiedNodeName     qnodename2 = qnodename.popFront();
    xmlpp::Node::NodeList children   = node->get_children();

    for(xmlpp::Node::NodeList::const_iterator child_itr = children.begin(); child_itr != children.end(); ++child_itr)
    {
      if((*child_itr)->get_name() == qnodename2.getFirstName())
      {
        try
        {
          return getNode(qnodename2, *child_itr);
        }
        catch(const std::exception & e)
        { }
      }
    }

    throw std::exception();
  }
}

//==========================================================
//
//  getNodes(...)
//
//==========================================================
inline XMLParser::NodePtrVector
XMLParser::getNodes(const QualifiedNodeName & qnodename, NodePtr _node) const
{
  // If the default argument (NULL) was used for _node, then use this object's root node:
  
  NodePtr node = _node == NULL ? getRootNode() : _node;

  // Set up return vector:
  
  NodePtrVector node_list;

  // Is the first name correct?

  if(qnodename.getNodeCount() == 1)
  {
    if(node->get_name() == qnodename.getFirstName())
      node_list.push_back(node);
  }
  else // nodecount > 1
  {
    QualifiedNodeName     qnodename2 = qnodename.popFront();
    xmlpp::Node::NodeList children   = node->get_children();

    for(xmlpp::Node::NodeList::const_iterator child_itr = children.begin(); child_itr != children.end(); ++child_itr)
    {
      NodePtrVector child_node_list = getNodes(qnodename2, *child_itr);
      
      for(size_t i = 0; i < child_node_list.size(); ++i)
        node_list.push_back(child_node_list[i]);
    }
  }
  
  return node_list;
}

//=========================================================
//
//  extractTextContent(...)
//
//=========================================================
inline std::string
XMLParser::extractTextContent(NodePtr node) const
{
  xmlpp::Node::NodeList children = node->get_children();

  for(xmlpp::Node::NodeList::const_iterator children_itr = children.begin(); children_itr != children.end(); ++children_itr)
  {
    const xmlpp::TextNode * node_text = dynamic_cast<const xmlpp::TextNode*> (*children_itr);
          
    if(node_text)
      return node_text->get_content();
  }

  return "";
}   
                          
inline XMLParser::AttributeList 
XMLParser::getAttributes(NodePtr node) const
{
  if(const xmlpp::Element* node_element = dynamic_cast<const xmlpp::Element*>(node))
  {
    return node_element->get_attributes();
  }
  else
  {
    return AttributeList();
  }
}

//================================================
//
//  getRootNode()
//
//================================================
inline XMLParser::NodePtr
XMLParser::getRootNode() const
{
  return my_parser.get_document()->get_root_node();
}

//===============================================
//
//  toString()
//
//===============================================
inline std::string 
XMLParser::toString() const
{
  // You have to const_cast the my_parser member because write_to_string() is non-const.
  
  return const_cast<xmlpp::DomParser*>(&my_parser)->get_document()->write_to_string();
}

} // End namespace ERROR
} // End namespace SAGE

#endif
