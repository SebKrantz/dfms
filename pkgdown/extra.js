// Process raw LaTeX math and render with KaTeX
// Handles both inline math \(...\) and display math $$...$$
(function() {
  function renderAllMath() {
    if (typeof katex === 'undefined') {
      setTimeout(renderAllMath, 50);
      return;
    }
    
    // Helper function to find matching closing delimiter
    function findMatchingClose(text, start, openDelim, closeDelim) {
      var depth = 1;
      var i = start + openDelim.length;
      while (i < text.length && depth > 0) {
        if (text.substring(i, i + openDelim.length) === openDelim) {
          depth++;
          i += openDelim.length;
        } else if (text.substring(i, i + closeDelim.length) === closeDelim) {
          depth--;
          if (depth === 0) {
            return i + closeDelim.length;
          }
          i += closeDelim.length;
        } else {
          i++;
        }
      }
      return -1;
    }
    
    // Process inline math: \(...\)
    function processInlineMath(node) {
      var text = node.textContent;
      var result = [];
      var i = 0;
      
      while (i < text.length) {
        var openIdx = text.indexOf('\\(', i);
        if (openIdx === -1) {
          // No more math, add remaining text
          if (i < text.length) {
            result.push({type: 'text', content: text.substring(i)});
          }
          break;
        }
        
        // Add text before the math
        if (openIdx > i) {
          result.push({type: 'text', content: text.substring(i, openIdx)});
        }
        
        // Find matching closing delimiter
        var closeIdx = findMatchingClose(text, openIdx, '\\(', '\\)');
        if (closeIdx === -1) {
          // No matching close, treat as regular text
          result.push({type: 'text', content: text.substring(openIdx)});
          break;
        }
        
        // Extract math content (without delimiters)
        var mathContent = text.substring(openIdx + 2, closeIdx - 2);
        result.push({type: 'math', content: mathContent, display: false});
        i = closeIdx;
      }
      
      return result;
    }
    
    // Process all text nodes
    var textNodes = [];
    var walker = document.createTreeWalker(
      document.body,
      NodeFilter.SHOW_TEXT,
      {
        acceptNode: function(node) {
          // Skip script and style nodes
          var parent = node.parentNode;
          if (parent && (parent.tagName === 'SCRIPT' || parent.tagName === 'STYLE')) {
            return NodeFilter.FILTER_REJECT;
          }
          return node.textContent.includes('\\(') ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_REJECT;
        }
      },
      false
    );
    
    var node;
    while (node = walker.nextNode()) {
      textNodes.push(node);
    }
    
    textNodes.forEach(function(textNode) {
      var parent = textNode.parentNode;
      if (!parent) return;
      
      var parts = processInlineMath(textNode);
      if (parts.length <= 1 && parts[0] && parts[0].type === 'text') {
        // No math found, skip
        return;
      }
      
      var fragment = document.createDocumentFragment();
      parts.forEach(function(part) {
        if (part.type === 'text') {
          fragment.appendChild(document.createTextNode(part.content));
        } else if (part.type === 'math') {
          var span = document.createElement('span');
          span.className = 'math inline';
          try {
            katex.render(part.content, span, {throwOnError: false, displayMode: part.display});
            fragment.appendChild(span);
          } catch (e) {
            // If rendering fails, show original
            fragment.appendChild(document.createTextNode('\(' + part.content + '\)'));
          }
        }
      });
      
      parent.replaceChild(fragment, textNode);
    });
    
    // Process display math: $$...$$
    var allElements = document.querySelectorAll('p, td, li, dd, dt');
    allElements.forEach(function(el) {
      var html = el.innerHTML;
      if (html.includes('$$')) {
        var newHTML = html.replace(/\$\$([\s\S]*?)\$\$/g, function(match, math) {
          var span = document.createElement('span');
          span.className = 'math display';
          try {
            katex.render(math.trim(), span, {throwOnError: false, displayMode: true});
            return span.outerHTML;
          } catch (e) {
            return match;
          }
        });
        if (newHTML !== html) {
          el.innerHTML = newHTML;
        }
      }
    });
  }
  
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', function() {
      setTimeout(renderAllMath, 200);
    });
  } else {
    setTimeout(renderAllMath, 200);
  }
})();
