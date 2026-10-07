// Panel behaviour: the interaction builder's pool, and one row per factor in
// Reference Levels and per linear covariate in Covariate Scaling.
// Adapted from the prognosis module's cox.events.js (jus 2.0 events API).

const events = {
    update: function(ui) {
        calcInteractionPool(ui, this);
        updateRefLevels(ui, this);
        updateLevelControls(ui, this);
        updateScales(ui, this);
    },

    onChange_predictors: function(ui) {
        calcInteractionPool(ui, this);
        updateRefLevels(ui, this);
        updateLevelControls(ui, this);
        updateScales(ui, this);
    },

    onChange_refLevels: function(ui) {
        updateLevelControls(ui, this);
    },

    onUpdate_interactionSupplier: function(ui) {
        let vars = collectPredictors(ui, this);
        ui.interactionSupplier.setValue(this.valuesToItems(vars, FormatDef.variable));
    }
};

let collectPredictors = function(ui, context) {
    let a = context.cloneArray(ui.covs.value(), []);
    let b = context.cloneArray(ui.splineVars.value(), []);
    let c = context.cloneArray(ui.factors.value(), []);
    return a.concat(b).concat(c);
};

// terms that reference a removed predictor are pruned
let calcInteractionPool = function(ui, context) {
    let vars = collectPredictors(ui, context);
    ui.interactionSupplier.setValue(context.valuesToItems(vars, FormatDef.variable));

    let varsDiff = context.findChanges("predictorList", vars, true, FormatDef.variable);
    let termsList = context.cloneArray(ui.interactions.value(), []);
    let changed = false;
    for (let i = 0; i < varsDiff.removed.length; i++) {
        for (let j = 0; j < termsList.length; j++) {
            if (FormatDef.term.contains(termsList[j], varsDiff.removed[i])) {
                termsList.splice(j, 1);
                changed = true;
                j -= 1;
            }
        }
    }
    if (changed)
        ui.interactions.setValue(termsList);
};

// keep the choice already made for a variable that stays in the list
let syncList = function(list, vars, make) {
    let out = [];
    for (let i = 0; i < vars.length; i++) {
        let found = null;
        for (let j = 0; j < list.length; j++) {
            if (list[j].var === vars[i]) {
                found = list[j];
                break;
            }
        }
        out.push(found === null ? make(vars[i]) : found);
    }
    return out;
};

let updateRefLevels = function(ui, context) {
    let factors = context.cloneArray(ui.factors.value(), []);
    let current = context.cloneArray(ui.refLevels.value(), []);
    ui.refLevels.setValue(syncList(current, factors, (v) => ({ var: v, ref: null })));
};

// each LevelSelector lists the levels of its own row's variable
let updateLevelControls = function(ui, context) {
    let list = ui.refLevels.value();
    ui.refLevels.applyToItems(0, (item, index, column) => {
        if (column === 1)
            item.setPropertyValue('variable', list[index].var);
    });
};

let updateScales = function(ui, context) {
    let covs = context.cloneArray(ui.covs.value(), []);
    let current = context.cloneArray(ui.covScales.value(), []);
    ui.covScales.setValue(syncList(current, covs, (v) => ({ var: v, scale: 'unit' })));
};

module.exports = events;
