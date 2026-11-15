{{ name | escape | underline}}


.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :exclude-members: _validate_grammar, CubicSpline{% if objname == 'Arrow3D' %}, set{% endif %}

{# Inherited members are linked below but must NOT be autodoc-registered on
   every subclass. They are documented on their defining class instead.
   Enabling inherited-member registration here creates duplicate descriptions
   for Repeated, MemberChange, PyParsing classes, etc. #}

{% if not objname.startswith('_') %}
   .. rubric:: {{ _('Class hierarchy') }}
   .. inheritance-diagram:: {{ fullname }}
      :parts: 2
{% endif %}

   |

   .. rubric:: {{ _('Constructor') }}
   {% block init %}
   .. automethod:: __init__
   {% endblock %}

   {% block own_attributes %}
   {% if own_attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in own_attributes if item not in ('_validate_grammar', 'CubicSpline') %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block inherited_attributes %}
   {% if inherited_attributes and objname != 'Arrow3D' %}
   .. rubric:: {{ _('Inherited attributes') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in inherited_attributes if item not in ('_validate_grammar', 'CubicSpline') %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block own_methods %}
   {% if own_methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in own_methods if item != '_validate_grammar' and (objname != 'Arrow3D' or item != 'set') %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block inherited_methods %}
   {% if inherited_methods and objname != 'Arrow3D' %}
   .. rubric:: {{ _('Inherited methods') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in inherited_methods if item != '_validate_grammar' and (objname != 'Arrow3D' or item != 'set') %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}
