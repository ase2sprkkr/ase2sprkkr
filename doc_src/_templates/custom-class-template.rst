{{ name | escape | underline}}


.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :inherited-members:

   .. rubric:: {{ _('Class hierarchy') }}
   .. inheritance-diagram:: {{ fullname }}
      :parts: 2

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
     {% for item in own_attributes %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block inherited_attributes %}
   {% if inherited_attributes %}
   .. rubric:: {{ _('Inherited attributes') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in inherited_attributes %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block own_methods %}
   {% if own_methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in own_methods %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block inherited_methods %}
   {% if inherited_methods %}
   .. rubric:: {{ _('Inherited methods') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in inherited_methods %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}
