{{ name | escape | underline}}

Full name: **{{ fullname | escape }}**

{# Sphinx excludes private classes from inheritance diagrams by default.
   Do not generate an empty diagram for a module with only private classes. #}
{% set hierarchy_classes = [] %}
{% for item in classes if not item.startswith('_') %}
{% set _ = hierarchy_classes.append(fullname ~ '.' ~ item) %}
{% endfor %}
{% if hierarchy_classes %}

.. rubric:: Module class hierarchy

.. inheritance-diagram:: {{ hierarchy_classes | join(' ') }}
      :parts: 1

|
{% endif %}

.. rubric:: Description

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module Attributes

   .. autosummary::
      :toctree:
      :template: custom-base-template.rst
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
      :template: custom-base-template.rst
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :template: custom-class-template.rst
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


  {% block modules %}
  {% if modules %}
  .. rubric:: Modules

  .. autosummary::
     :toctree:
     :template: custom-module-template.rst
     :recursive:
  {% for item in modules %}
  {% if item != 'test' and item != 'tests' and '.test.' not in item and '.tests.' not in item and not item.endswith('.test') and not item.endswith('.tests') %}
     {{ item }}
  {% endif %}
  {%- endfor %}
  {% endif %}
  {% endblock %}
